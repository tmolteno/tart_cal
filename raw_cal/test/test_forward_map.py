#
# # Build a sky model with some sources at known locations (el, az)
# # Build a telescope model with some positions (x,y)
# # Construct an image of the sky using inverse FFT
# # Find the peaks, and convert the image coordinates to el,az
# # Assert that the el, az match.
# #
# #
# import numpy as np
# import torch
#
# class Source:
#     def __init__(self, el, az, jy):
#         self.el = el
#         self.az = az
#         self.jy = jy
#
# class Telescope:
#
#     def __init__(self, ant_pos):
#         self.ant_pos = ant_pos
#         self.nant = len(ant_pos)
#         self.baseline_indices =
#
#         baselines = self.get_baselines(ant_pos)
#         self.max_bl = np.max(np.abs(baselines))
#
#         num_bin = 2**9  # Image resolution
#         nw=num_bin/4
#
#         self.uv_grid = torch.zeros((num_bin, num_bin))
#
#     def get_baselines(self, ant_pos):
#         self.baselines = [(ant_pos[i] - ant_pos[j]) for i,j in self.baseline_indices]
#
#     def get_vis(self, sources):
#         # form UVW grid
#         x = torch.rand(10, 10, dtype=torch.complex64)
#         fft2 = torch.fft.fft2(x)
#
#     def get_image_from_vis(vis, bl):
#
#         uv_plane = self.uv_grid * 0.0
#         for v in zip(vis, bl):
#             uv_plane[]
#         cal_ift = fft.fftshift(fft.ifft2(fft.ifftshift(uv_grid)))
#
#         # Take the absolute value to make an intensity image
#         img = np.abs(cal_ift)
#
# def do_image(uv_grid):
#     cal_ift = fft.fftshift(fft.ifft2(fft.ifftshift(uv_grid)))
#     img = torch.abs(cal_ift)
#
# def forward_map(uv_grid, baselines, vis, rotation_angle, complex_gains):
#     vis = vis * complex_gains   # Elementwise product of gains with visibilities
#     rotated_pos = torch.rotate(ant_pos, rotation_angle)
#     return do_image(uv_grid)
