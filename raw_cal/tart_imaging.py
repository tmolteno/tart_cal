import numpy as np
import torch

from tart.imaging import elaz

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
        self.pixels[x, y] += power

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


def get_uv_array_indices(image_size, baselines, wavelength):
    ''' A little function to produce the index into the u-v array
        for a given value (u, measured in wavelengths)
    '''
    uv_max = image_size / 4   # (1.2 * np.pi)
    middle = image_size / 2

    uv = baselines[:, 0:2] / wavelength

    uv_pix = middle + (uv / uv_max)*(image_size/2)

    # Truncate any uv elements that don't fit into the array
    max_pix = np.max(uv_pix, axis=1)
    uv_pix = uv_pix[max_pix < image_size, :]

    return np.round(uv_pix).astype(int)

def get_simulated_visibilities(image, baselines, wavelength):
    ''' Simulate some visibilities from an image by
    doing an fft and then sampling at the correct points
    and generate a gridded plane'''

    uv_image = get_uv_plane(image)
    imsize = uv_image.image_size

    uvpix = get_uv_array_indices(imsize, baselines, wavelength=wavelength)
    upix = uvpix[:,0]
    vpix = uvpix[:,1]
    vis = uv_image.pixels[upix, vpix]

    return vis


def gridder(vis, imsize, baselines, wavelength):
    gridded = Image.zeros(imsize)

    uvpix = get_uv_array_indices(imsize, baselines, wavelength=wavelength)
    upix = uvpix[:,0]
    vpix = uvpix[:,1]
    gridded.pixels[upix,vpix] += vis

    uvpix_neg = get_uv_array_indices(imsize, -baselines, wavelength=wavelength)
    upix = uvpix_neg[:,0]
    vpix = uvpix_neg[:,1]
    gridded.pixels[upix,vpix] += torch.conj(vis)

    return gridded
