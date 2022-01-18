/* 
 * Utility to convert photos into a Ghibli style image.
 *
 * 
 * Copyright (C) 2022 Peter Pham
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), 
 * to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, 
 * and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
 * 
 */

#include <stdlib.h>
#include <stdio.h>
#include <memory.h>
#include <iostream>
#include <string>
#include <vector>
#include <png.h>
#include <math.h> 

using namespace std;


void throw_error(const char* estr)
{
    cout << estr << endl;
    abort();
}

int image_width;
int image_height;

png_infop info_ptr;
png_bytepp idata;
int bpp;

void read_png(const char* filename) 
{
    /* Check if png */
    FILE *fptr = fopen(filename, "rb");
    if (!fptr) {
        throw_error("File could not be opened.");
    }

    char header[8];
    fread(header, 1, 8, fptr);
    if (png_sig_cmp((png_const_bytep)header, 0, 8)) {
        throw_error("File is not recognized as a png.");
    }

    /* Set up and read */
    png_structp png_ptr = png_create_read_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
    info_ptr = png_create_info_struct(png_ptr);
    png_init_io(png_ptr, fptr);
    png_set_sig_bytes(png_ptr, 8);

    png_read_png(png_ptr, info_ptr, PNG_TRANSFORM_IDENTITY, NULL);
    idata = png_get_rows(png_ptr, info_ptr);

    image_width = png_get_image_width(png_ptr, info_ptr);
    image_height = png_get_image_height(png_ptr, info_ptr);
    bpp = (png_get_color_type(png_ptr, info_ptr) == PNG_COLOR_TYPE_RGBA) ? 4 : 3;

    png_destroy_read_struct(&png_ptr, NULL, NULL);
    fclose(fptr);
}

void write_png(const char* filename, png_bytepp data)
{
    /* Check if png */
    FILE *fptr = fopen(filename, "wb");
    if (!fptr) {
        throw_error("File could not be opened/created.");
    }

    /* Set up and write */
    png_structp png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
    png_init_io(png_ptr, fptr);
    png_set_rows(png_ptr, info_ptr, data);
    png_write_png(png_ptr, info_ptr, PNG_TRANSFORM_IDENTITY, NULL);

    png_destroy_write_struct(&png_ptr, &info_ptr);
    fclose(fptr);
}

int bound(int n, int low, int high) 
{
    return min(max(n, low), high);
}

float fbound(float n, float low, float high) 
{
    return min(max(n, low), high);
}

png_bytepp filter(png_bytepp img, vector<vector<double>> filter_m, const double factor) 
{
    // This is slow but yolo (only dealing with small convolution matrix anyways)
    png_bytepp out = (png_bytepp)malloc(sizeof(png_bytep)*image_height);
    for(int i = 0; i < image_height; i++) {
        out[i] = (png_bytep)malloc(image_width*bpp);
    }

    int filter_width = filter_m[0].size();
    int filter_height = filter_m.size();

    for(int y = 0; y < image_height; y++) {
        for(int x = 0; x < image_width; x++) 
        {
            double red = 0.0, green = 0.0, blue = 0.0;
        
            for(int f_y = 0; f_y < filter_height; f_y++) {
                for(int f_x = 0; f_x < filter_width; f_x++) 
                {
                    int fpixel_x = bound(x-filter_width/2+f_x, 0, image_width-1);
                    int fpixel_y = bound(y-filter_height/2+f_y, 0, image_height-1);

                    png_bytep px = &img[fpixel_y][fpixel_x*bpp];
                    red += px[0] * filter_m[f_y][f_x];
                    green += px[1] * filter_m[f_y][f_x];
                    blue += px[2] * filter_m[f_y][f_x];
                }
            }

            png_bytep outpx = &out[y][x*bpp];
            outpx[0] = bound(int(red * factor), 0, 255);
            outpx[1] = bound(int(green * factor), 0, 255);
            outpx[2] = bound(int(blue * factor), 0, 255);
            if (bpp == 4)
                outpx[3] = (&img[y][x*bpp])[3];
        }
    }

    return out;
}

png_bytepp gaussian_b(png_bytepp img, const int amnt) 
{
    vector<vector<double>> filter_m;
    // {
    //     {0, 0, 1, 0, 0},
    //     {0, 1, 2, 1, 0},
    //     {1, 2, 3, 2, 1},
    //     {0, 1, 2, 1, 0},
    //     {0, 0, 1, 0, 0},
    // }; etc.

    for(int i = 0; i < amnt; i++) {
        vector<double> row;
        for(int j = 0; j < amnt; j++) {
            row.push_back(max(0, (amnt/2+1)-(abs(i-amnt/2)+abs(j-amnt/2))));
        }
        filter_m.push_back(row);
    }

    const double factor = 1.0 / double((amnt/2+1)*(2*(amnt/2+1)*(amnt/2+1)+1)/3); 

    return filter(img, filter_m, factor);
}

/* Compute luminance of rgb values according to D50 spec. */
double luminance(int r, int g, int b)
{
    const double RGB_LUM_R = 0.22245;
    const double RGB_LUM_G = 0.71695;
    const double RGB_LUM_B = 0.06070;

    return (r*RGB_LUM_R) + (g*RGB_LUM_G) + (b*RGB_LUM_B);
}

png_bytepp oilify(png_bytepp img, const int radius, const double intensity_lvls, const int exponent)
{
    /* Generate intensity map */
    vector<vector<int>> inten;

    for(int y = 0; y < image_height; y++) {
        vector<int> row;
        for(int x = 0; x < image_width; x++) 
        {
            png_bytep px = &img[y][x*bpp];
            int lum = (luminance(px[0], px[1], px[2]) * intensity_lvls) / 255;
            lum = bound(lum, 0, 255);
            row.push_back(lum);
        }
        inten.push_back(row);
    }

    /* Apply oilify effect */
    png_bytepp out = (png_bytepp)malloc(sizeof(png_bytep)*image_height);
    for(int i = 0; i < image_height; i++) {
        out[i] = (png_bytep)malloc(image_width*bpp);
    }

    for(int y = 0; y < image_height; y++) {
        for(int x = 0; x < image_width; x++) 
        {
            vector<int> hist_inten(256, 0);
            vector<vector<int>> hist_rgb(256, vector<int>{0, 0, 0, 255});

            for(int r_y = -radius; r_y <= radius; r_y++) {
                for(int r_x = -radius; r_x <= radius; r_x++) 
                {
                    // circular mask
                    if (r_y*r_y + r_x*r_x > radius*radius) 
                        continue;
                    int i_y = bound(y+r_y, 0, image_height-1);
                    int i_x = bound(x+r_x, 0, image_width-1);

                    int cur_intensity = inten[i_y][i_x];
                    hist_inten[cur_intensity]++;

                    png_bytep px = &img[i_y][i_x*bpp];
                    hist_rgb[cur_intensity][0] += px[0];
                    hist_rgb[cur_intensity][1] += px[1];
                    hist_rgb[cur_intensity][2] += px[2];
                }
            }

            /* TODO: Weighted average color w exponent */
            int max_cnt = 0;
            int max_inten = 0;
            for(int i = 0; i <= 255; i++) {
                if (hist_inten[i] > max_cnt) {
                    max_cnt = hist_inten[i];
                    max_inten = i;
                }
            }

            png_bytep outpx = &out[y][x*bpp];
            outpx[0] = hist_rgb[max_inten][0] / max_cnt;
            outpx[1] = hist_rgb[max_inten][1] / max_cnt;
            outpx[2] = hist_rgb[max_inten][2] / max_cnt;
            if (bpp == 4)
                outpx[3] = (&img[y][x*bpp])[3];
        }
    }

    return out;
}

void rgb2hsv(png_bytep px, float* hsv)
{
    float rc = (float)px[0]/255.0f;
    float gc = (float)px[1]/255.0f;
    float bc = (float)px[2]/255.0f;

    float cmax = max(max(rc, gc), bc);
    float cmin = min(min(rc, gc), bc);
    float diff = cmax-cmin;

    if (diff == 0)
        hsv[0] = 0.0f;
    else if (cmax == rc)
        hsv[0] = fmodf((60.0f*((gc-bc)/diff)+360.0f), 360.0f);
    else if (cmax == gc)
        hsv[0] = fmodf((60.0f*((bc-rc)/diff)+120.0f), 360.0f);
    else if (cmax == bc)
        hsv[0] = fmodf((60.0f*((rc-gc)/diff)+240.0f), 360.0f);

    if (cmax == 0)
        hsv[1] = 0;
    else
        hsv[1] = (diff/cmax)*100.0f;

    hsv[2] = cmax*100.0f;
}

void hsv2rgb(float* hsv, png_bytep px)
{
    float ss = hsv[1]/100.0f;
    float vv = hsv[2]/100.0f;

    if (hsv[0] >= 360.0f)
        hsv[0] = 0.0f;
    hsv[0] /= 60.0f;
    long i(hsv[0]);
    float ff = hsv[0]-i;

    float p = vv*(1.0f-ss);
    float q = vv*(1.0f-(ss*ff));
    float t = vv*(1.0f-(ss*(1.0f-ff)));

    switch(i) 
    {
        case 0:
            px[0] = vv*255;
            px[1] = t*255;
            px[2] = p*255;
            break;
        case 1:
            px[0] = q*255;
            px[1] = vv*255;
            px[2] = p*255;
            break;
        case 2:
            px[0] = p*255;
            px[1] = vv*255;
            px[2] = t*255;
            break;
        case 3:
            px[0] = p*255;
            px[1] = q*255;
            px[2] = vv*255;
            break;
        case 4:
            px[0] = t*255;
            px[1] = p*255;
            px[2] = vv*255;
            break;
        default:
            px[0] = vv*255;
            px[1] = p*255;
            px[2] = q*255;
            break;
    }
}

png_bytepp hsvcorrect_noise(png_bytepp img)
{
    /* Convert to hsv */
    float** img_hsv = (float**)malloc(sizeof(float)*image_width*image_height*3);
    for(int i = 0; i < image_height; i++) {
        img_hsv[i] = (float*)malloc(image_width*3*sizeof(float));
    }

    float avg_hsv[3] = {0.0f, 0.0f, 0.0f};
    for(int y = 0; y < image_height; y++) {
        vector<float*> row;
        for(int x = 0; x < image_width; x++) 
        {
            float* hsv = &img_hsv[y][x*3];
            rgb2hsv(&img[y][x*bpp], hsv);

            avg_hsv[0] += hsv[0];
            avg_hsv[1] += hsv[1];
            avg_hsv[2] += hsv[2];
            if (x == 0 && y ==0)
                cout << hsv[0] << "," << hsv[1] << "," << hsv[2] << endl;
        }
    }
    float px_cnt = image_height*image_width;
    avg_hsv[0] /= px_cnt;
    avg_hsv[1] /= px_cnt;
    avg_hsv[2] /= px_cnt;
    cout << avg_hsv[0] << "," << avg_hsv[1] << "," << avg_hsv[2] << endl;

    /* Determine changes */
    float hue_ratio = 1.0f, sat_ratio = 1.0f, val_ratio = 1.0f;

    if (avg_hsv[0] < 340 && avg_hsv[0] > 60 && avg_hsv[2] < 76) {
        if (avg_hsv[1] < 90) {
            sat_ratio = 1.0f-(90.0f-avg_hsv[1])/360.0f;
            val_ratio = 1.0f-(80.0f-avg_hsv[1])/640.0f;
        }
        if (avg_hsv[0] < 216) {
            hue_ratio = 1.0f+(216.0f-avg_hsv[0])/720.0f;
        } else {
            hue_ratio = 1.0f-(avg_hsv[0]-216.0f)/960.0f;
        }
    } else {
        if (avg_hsv[1] < 80) {
            sat_ratio = 1.0f+(100.0f-avg_hsv[1])/500.0f;
        }
        if (avg_hsv[0] < 32) {
            hue_ratio = 1.0f-(avg_hsv[0])/96.0f;
        } else if (avg_hsv[0] > 254) {
            hue_ratio = 1.0f+(360.0f-avg_hsv[0])/1080.0f;
        }
    }
    cout << avg_hsv[0] << "," << avg_hsv[1] << "," << avg_hsv[2] << endl;
    cout << hue_ratio << "," << sat_ratio << "," << val_ratio << endl;
    cout << (&img_hsv[100][300])[0] << "," << (&img_hsv[100][300])[1] << "," << (&img_hsv[100][300])[2] << endl;

    /* Apply changes and noise */
    int hue_noise = 5;
    int sat_noise = 4;
    int val_noise = 9;

    png_bytepp out = (png_bytepp)malloc(sizeof(png_bytep)*image_height);
    for(int i = 0; i < image_height; i++) {
        out[i] = (png_bytep)malloc(image_width*bpp);
    }

    for(int y = 0; y < image_height; y++) {
        for(int x = 0; x < image_width; x++) 
        {
            float* hsvpx = (&img_hsv[y][x*3]);
            hsvpx[0] = fbound(hsvpx[0]*hue_ratio+(rand()%(hue_noise*2)-hue_noise), 0.0f, 360.0f);
            hsvpx[1] = fbound(hsvpx[1]*sat_ratio+(rand()%(sat_noise*2)-sat_noise), 0.0f, 100.0f);
            hsvpx[2] = fbound(hsvpx[2]*val_ratio+(rand()%(val_noise*2)-val_noise), 0.0f, 100.0f);
            /* HSV noise */
            

            png_bytep outpx = &out[y][x*bpp];
            hsv2rgb(hsvpx, outpx);
        }
    }

    free(img_hsv);

    return out;
}

png_bytepp exposure(png_bytepp img, const double exposure, const double black)
{
    png_bytepp out = (png_bytepp)malloc(sizeof(png_bytep)*image_height);
    for(int i = 0; i < image_height; i++) {
        out[i] = (png_bytep)malloc(image_width*bpp);
    }

    for(int y = 0; y < image_height; y++) {
        for(int x = 0; x < image_width; x++) 
        {
            double expos = exp2(-exposure);
            double gain = 1.0/(max(expos-black, 0.001));

            png_bytep px = &img[y][x*bpp];
            png_bytep outpx = &out[y][x*bpp];
            outpx[0] = bound((int)(((double)px[0]-black*255.0)*gain), 0, 255);
            outpx[1] = bound((int)(((double)px[1]-black*255.0)*gain), 0, 255);
            outpx[2] = bound((int)(((double)px[2]-black*255.0)*gain), 0, 255);
            if (bpp == 4)
                outpx[3] = (&img[y][x*bpp])[3];
        }
    }

    return out;
}

png_bytepp ghiblify(png_bytepp img) 
{
    png_bytepp out;

    const int blur_factor1 = max((image_width+image_height)/610 + 1, 3);
    const int blur_factor2 = max((image_width+image_height)/1024 + 1, 2);
    const int oil_radius = max((image_width+image_height)/900 + 1, 2);
    const double exposure_level = 0.4;
    const double black_level = 0.05;

    out = gaussian_b(img, blur_factor1);
    out = oilify(out, oil_radius, 64, 2);
    out = gaussian_b(out, blur_factor2);
    out = exposure(out, exposure_level, black_level);
    out = hsvcorrect_noise(out);

    return out;
}

int main(int argc, char *argv[])
{
    float hsv[3];

    // unsigned char b[4] = {124u,168u,71u,255u}; 
    // cout << (int)b[0] << "," <<  (int)b[1] << "," << (int)b[2] << endl;
    // rgb2hsv((png_bytep)b, hsv);
    // cout << hsv[0] << "," <<  hsv[1] << "," << hsv[2] << endl;
    // hsv2rgb(hsv, (png_bytep)b);
    // cout << (int)b[0] << "," <<  (int)b[1] << "," << (int)b[2] << endl;
    read_png(argv[1]);

    string out;
    if (argc < 3) {
        string fn = string(argv[1]);
        out = fn.substr(0, fn.length()-4)+"_ghibli.png";
    } else {
        out = argv[2];
    }
    write_png(out.c_str(), ghiblify(idata));
    cout << "Ghiblified image outputted to " << out << endl;
    
    return 0;
}