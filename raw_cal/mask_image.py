import logging

from tart.imaging import imaging

logger = logging.getLogger()


def add_source(mask, src, radius_deg):
    '''
        mask: array-like
        src: l, m coordinates
    '''
    image_size = mask.shape[0]

    logger.info(f"    add_source {src}")

    x0, y0 = src.get_px(image_size)
    r_pix = imaging.deg_to_pix(image_size, radius_deg)**2
    for y in range(image_size):
        for x in range(image_size):
            r2 = (y - y0)**2 + (x - x0)**2
            if r2 <= r_pix:
                mask[x, y] = 1.0

