import unittest
import torch
import matplotlib.pyplot as plt
import numpy as np

from raw_cal import tart_imaging


class TestImaging(unittest.TestCase):

    def setUp(self):
        self.tart_wavelength = 2.9979e8 / 1.57542e9

    def random_image(self, n_points, image_size):
        larray = np.random.uniform(-1, 1, size=n_points)
        marray = np.random.uniform(-1, 1, size=n_points)
        parray = np.random.uniform(5, 10, size=n_points)

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

    def test_simulated_vis(self):
        imsize = 128
        npoint = 10
        img, larray, marray, parray = self.random_image(npoint, image_size=imsize)

        # Now generate a random antenna array
        n_ant = 24
        ant_pos = self.random_ant_pos(n_ant)
        bl = tart_imaging.get_baselines(ant_pos)

        vis = tart_imaging.get_simulated_visibilities(
            image=img,
            baselines=bl, wavelength=self.tart_wavelength)

        uv_image = tart_imaging.get_uv_plane(img)
        uvpix = tart_imaging.get_uv_array_indices(imsize, bl, wavelength=self.tart_wavelength)

        for i in range(len(vis)):
            u, v = uvpix[i]
            self.assertAlmostEqual(vis[i],  uv_image.pixels[u, v])

    def test_visibility_conjugate(self):
        imsize = 128
        npoint = 10
        img, larray, marray, parray = self.random_image(npoint, image_size=imsize)

        # Now generate a random antenna array
        n_ant = 24
        ant_pos = self.random_ant_pos(n_ant)
        bl = tart_imaging.get_baselines(ant_pos)

        wavelength = self.tart_wavelength

        vis = tart_imaging.get_simulated_visibilities(img, bl, wavelength=wavelength)
        uvpix = tart_imaging.get_uv_array_indices(imsize, bl, wavelength=wavelength)

        vis_neg = tart_imaging.get_simulated_visibilities(img, -bl, wavelength=wavelength)
        uvpix_neg = tart_imaging.get_uv_array_indices(imsize, -bl, wavelength=wavelength)

        for vp, vn, uv, uvn in zip(vis, vis_neg, uvpix, uvpix_neg):
            print(vp, vn, uv, uvn)
            self.assertAlmostEqual(vp.item(), torch.conj(vn).item(), 4)

    def test_multipoint_imaging(self):
        imsize = 256
        npoint = 10
        img, larray, marray, parray = self.random_image(npoint, image_size=imsize)

        # Now generate a random antenna array
        n_ant = 24
        ant_pos = self.random_ant_pos(n_ant)
        bl = tart_imaging.get_baselines(ant_pos)

        # Simulate some visibilities
        vis = tart_imaging.get_simulated_visibilities(image=img, baselines=bl, wavelength=self.tart_wavelength)
        gridded = tart_imaging.gridder(vis=vis, image_size=imsize, baselines=bl, wavelength=self.tart_wavelength)

        uv_image = tart_imaging.get_uv_plane(img)

        # Image the masked array
        im2 = tart_imaging.get_image(gridded)

        # Scale the image as we haven't all visibilities
        tot_uv = torch.sum(torch.abs(uv_image.pixels))
        tot_masked = torch.sum(torch.abs(gridded.pixels))
        im2.pixels *= tot_uv / tot_masked

        if False:
            random_pix = np.abs(img.pixels.numpy())
            uv_pix = np.abs(uv_image.pixels.numpy())
            masked_pix = np.abs(gridded.pixels.numpy())
            recovered = np.abs(im2.pixels.numpy())

            fig, ax = plt.subplots(nrows=2, ncols=2)
            im00 = ax[0, 0].imshow(random_pix)
            fig.colorbar(im00, ax=ax[0, 0])
            im01 = ax[0, 1].imshow(uv_pix)
            fig.colorbar(im01, ax=ax[0, 1])

            im10 = ax[1, 0].imshow(masked_pix)
            fig.colorbar(im10, ax=ax[1, 0])

            im11 = ax[1, 1].imshow(recovered)
            fig.colorbar(im11, ax=ax[1, 1])
            plt.show()

        for i in range(npoint):
            p2 = im2.get_point(larray[i], marray[i])
            p2_scaled = torch.abs(p2).item()
            self.assertAlmostEqual(p2_scaled/parray[i], 1, delta=0.3)
