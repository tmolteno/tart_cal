import unittest

import matplotlib.pyplot as plt
import numpy as np

from raw_cal import mask_image
from tart.imaging import elaz

class TestMask(unittest.TestCase):

    def test_mask_diameter(self):

        image_pix = 128
        mask = np.zeros((image_pix, image_pix), dtype=np.float32)

        src = elaz.ElAz(90, 0)
        mask_image.add_source(mask, src, 2)

        self.assertAlmostEqual(mask[0,0], 0)

        center = image_pix // 2
        self.assertAlmostEqual(mask[center, center], 1.0)

        # plt.imshow(mask)
        # plt.show()
        #
        # raise RuntimeError("EEK")


    def test_mask_full(self):

        image_pix = 256
        mask = np.zeros((image_pix, image_pix), dtype=np.float32)

        src = elaz.ElAz(90, 0)
        mask_image.add_source(mask, src, 90)

        center = image_pix // 2
        self.assertAlmostEqual(mask[center, center], 1.0)
        #
        # plt.imshow(mask)
        # plt.show()
        #
        # raise RuntimeError("EEK")
