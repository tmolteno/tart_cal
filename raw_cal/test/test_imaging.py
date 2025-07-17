import unittest
import torch
import matplotlib.pyplot as plt
import numpy as np

from raw_cal import tart_imaging


class TestImaging(unittest.TestCase):

    def random_image(self, n_points, image_size):
        larray = np.random.uniform(-1, 1, size=n_points)
        marray = np.random.uniform(-1, 1, size=n_points)
        parray = np.random.uniform(1, 10, size=n_points)

        img = tart_imaging.Image.zeros(image_size)

        for i in range(n_points):
            img.add_point(larray[i], marray[i], power=parray[i])

        return img, larray, marray, parray

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

    def test_multipoint_imaging(self):
        tart_wavelength = 2.9979e8 / 1.57542e9
        imsize = 128
        npoint = 10
        img, larray, marray, parray = self.random_image(npoint, image_size=imsize)

        # Get the UV plane
        uv = tart_imaging.get_uv_plane(img)

        # Now generate a random antenna array
        n_ant = 24
        ant_pos = self.random_ant_pos(n_ant)

        # Now mask the uvplane with the antenna positions
        uvmasked = tart_imaging.Image.zeros(imsize)

        bl = tart_imaging.get_baselines(ant_pos)
        uvpix = tart_imaging.get_uv_coordinates(bl, imsize, wavelength=tart_wavelength)

        for u, v in uvpix:
            uvmasked.pixels[u, v] = uv.pixels[u, v]

        uvpix = tart_imaging.get_uv_coordinates(-bl, imsize, wavelength=tart_wavelength)

        for u, v in uvpix:
            uvmasked.pixels[u, v] = uv.pixels[u, v]

        # Image the masked array
        im2 = tart_imaging.get_image(uvmasked)

        tot_uv = torch.abs(torch.sum(uv.pixels))
        tot_masked = torch.abs(torch.sum(uvmasked.pixels))

        output_scale = tot_uv / tot_masked

        if True:
            random_pix = np.abs(img.pixels.numpy())
            random_pix = random_pix

            uv_pix = np.abs(uv.pixels.numpy())
            uv_pix = uv_pix

            tot_uv = np.sum(uv_pix)

            masked_pix = np.abs(uvmasked.pixels.numpy())
            masked_pix = masked_pix

            tot_masked = np.sum(masked_pix)

            recovered = np.abs(im2.pixels.numpy())
            recovered = recovered * (tot_uv / tot_masked)

            fig, ax = plt.subplots(nrows=2, ncols=2)
            im00 = ax[0,0].imshow(random_pix)
            fig.colorbar(im00, ax=ax[0,0])
            im01 = ax[0,1].imshow(uv_pix)
            fig.colorbar(im01, ax=ax[0,1])

            im10 = ax[1,0].imshow(masked_pix)
            fig.colorbar(im10, ax=ax[1,0])

            im11 = ax[1,1].imshow(recovered)
            fig.colorbar(im11, ax=ax[1,1])
            plt.show()

        for i in range(npoint):
            p2 = im2.get_point(larray[i], marray[i])
            p2_scaled = torch.abs(p2)*output_scale
            self.assertAlmostEqual(p2_scaled.item(), parray[i], 1)
