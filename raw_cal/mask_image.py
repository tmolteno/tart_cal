
import numpy as np
from tart.imaging import imaging

def add_source(mask, src, radius_deg):
    '''
        mask: array-like
        src: l, m coordinates
    '''
    image_size = mask.shape[0]

    thresh = 0.01
    y0, x0 = src.get_px(image_size)
    r_pix = imaging.deg_to_pix(image_size, radius_deg)**2
    for y in range(image_size):
        for x in range(image_size):
            r2 = (y - y0)**2 + (x - x0)**2
            if r2 > r_pix:
                p = 0
            else:
                p = 1.0
            mask[x, y] = p

