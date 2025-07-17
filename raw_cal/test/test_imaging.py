import unittest

import matplotlib.pyplot as plt
import numpy as np

from raw_cal import tart_imaging


class TestImaging(unittest.TestCase):

    def random_image(self, n_points):
        larray = np.random.uniform(-1, 1, size=100)
        marray = np.random.uniform(-1, 1, size=100)

        img = tart_imaging.Image.zeros(256)

        for i in range(100):
            img.add_point(larray[i], marray[i], power=10)

        return img

    def random_ant_pos(self, n_ant):
        ant_pos = np.random.uniform(-2, 2, size=(n_ant, 3))
        ant_pos[:, 2] = 0.0
        return ant_pos

    def test_to_uv_and_back(self):

        for i in range(100):
            l = np.random.uniform(-1, 1)
            m = np.random.uniform(-1, 1)
            power = np.random.uniform(1, 10)

            im = tart_imaging.Image.zeros(256)
            im.add_point(l, m, power=power)

            uv = tart_imaging.get_uv_plane(im)

            im2 = tart_imaging.get_image(uv)

            self.assertAlmostEqual(im2.get_point(l, m).item(), power, 5)

    def test_baselines(self):
        n_ant = 24
        ant_pos = self.random_ant_pos(n_ant)

        baselines = tart_imaging.get_baselines(ant_pos)

        self.assertEqual(baselines.shape[0], n_ant*(n_ant-1) / 2)
        self.assertEqual(baselines.shape[1], 3)

    def test_random_up(self):

        n_points = 10
        larray = np.random.uniform(-1, 1, size=n_points)
        marray = np.random.uniform(-1, 1, size=n_points)
        powers = np.random.uniform(1, 10, size=n_points)

        imsize = 128
        img = tart_imaging.Image.zeros(imsize)

        for i in range(n_points):
            img.add_point(larray[i], marray[i], power=powers[i])

        random_pix = np.abs(img.pixels.numpy())
        random_pix = random_pix / np.max(random_pix)

        uv = tart_imaging.get_uv_plane(img)
        uv_pix = np.abs(uv.pixels.numpy())
        uv_pix = uv_pix / np.max(uv_pix)

        # Now mask the uvplane with the antenna positions
        n_ant = 24
        ant_pos = self.random_ant_pos(n_ant)

        uvmasked = tart_imaging.Image.zeros(imsize)

        bl = tart_imaging.get_baselines(ant_pos)
        uvpix = tart_imaging.get_uv_coordinate(bl, imsize, wavelength=0.21)

        for u, v in uvpix:
            uvmasked.pixels[u, v] = uv.pixels[u, v]

        uvpix = tart_imaging.get_uv_coordinate(-bl, imsize, wavelength=0.21)

        for u, v in uvpix:
            uvmasked.pixels[u, v] = uv.pixels[u, v]

        masked_pix = np.abs(uvmasked.pixels.numpy())
        masked_pix = masked_pix / np.max(masked_pix)

        im2 = tart_imaging.get_image(uvmasked)
        recovered = np.abs(im2.pixels.numpy())
        recovered = recovered / np.max(recovered)

        fig, ax = plt.subplots(nrows=2, ncols=2)
        im = ax[0,0].imshow(random_pix)
        cbar = plt.colorbar(im)
        im = ax[0,1].imshow(uv_pix)
        cbar = plt.colorbar(im)

        im = ax[1,0].imshow(masked_pix)
        cbar = plt.colorbar(im)

        im = ax[1,1].imshow(recovered)
        cbar = plt.colorbar(im)
        plt.show()

        raise RuntimeError("ack")




