import numpy as np
import torch


class Image:

    def __init__(self, pixels):
        self.image_size = pixels.shape[0]
        self.width = self.image_size // 2
        self.pixels = pixels

    @classmethod
    def zeros(cls, image_size):
        pix = torch.zeros(image_size, image_size, dtype=torch.complex64)
        return Image(pix)

    def get_l_index(self, lm):
        return int(self.width*(1 + lm))

    def get_m_index(self, lm):
        return int(self.width*(1 - lm))

    def add_point(self, l, m, power):
        x = self.get_l_index(l)
        y = self.get_m_index(m)
        self.pixels[x, y] = power

    def get_point(self, l, m):
        x = self.get_l_index(l)
        y = self.get_m_index(m)
        return self.pixels[x, y]


def get_uv_plane(image):
    uv = torch.fft.ifftshift(
            torch.fft.fft2(
                torch.fft.fftshift(image.pixels)
                )
            )
    return Image(uv)


def get_image(uvplane):
    pix = torch.fft.fftshift(
            torch.fft.ifft2(
                torch.fft.ifftshift(uvplane.pixels)
                )
            )
    return Image(pix)


def get_baseline_indices(num_ant):
    bl_indices = []
    for i in range(num_ant-1):
        for j in range(i + 1, num_ant):
            bl_indices.append([i, j])
    return np.array(bl_indices)


def get_baselines(ant_pos):
    num_ant = ant_pos.shape[0]
    bl_indices = get_baseline_indices(num_ant)
    i_indices = bl_indices[:, 0]
    j_indices = bl_indices[:, 1]

    return ant_pos[j_indices, :] - ant_pos[i_indices, :]


def get_uv_coordinates(bl, image_size, wavelength):
    ''' A little function to produce the index into the u-v array
        for a given value (u, measured in wavelengths)
    '''
    uv_max = image_size / 4   # (1.2 * np.pi)
    middle = image_size / 2

    uv = bl[:, 0:2] / wavelength
    print(uv)
    print(f"uvmax {uv_max}")
    print(middle)

    uv_pix = middle + (uv / uv_max)*(image_size/2)

    max_pix = np.max(uv_pix, axis=1)
    print(max_pix)
    # goodu = (np.abs(u_pix) < uv_max).all()
    # goodv = (np.abs(v_pix) < uv_max).all()
    # good_uv = np.logical_and(goodu, goodv)
    uv_pix = uv_pix[max_pix < image_size, :]

    return np.floor(uv_pix).astype(int)
