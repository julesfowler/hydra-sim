
from astropy.io import fits
import glob
import numpy as np
import matplotlib.pyplot as plt

from PIL import Image



def apply_binning(data, resolution, modes):

    # taken almost exactly from Maaike's notebook

    if resolution==modes:
        binned_data = data

    else:
        binned_data = np.zeros((modes, modes))

        binning_factor = int(resolution/modes)
        #print(binning_factor)
        #print(binning_factor)
        #print(resolution-binning_factor)
        #print(resolution, modes)
        #print(binning_factor, resolution+1, binning_factor)
        for i_index in np.arange(binning_factor, resolution+1, binning_factor):
            for j_index in np.arange(binning_factor, resolution+1, binning_factor):
                binned_col = np.mean(
                             data[i_index-binning_factor:i_index,
                                  j_index-binning_factor:j_index]
                             )
                #print(i_index, j_index)
                #print(int((i_index+1)/binning_factor-1), int((j_index+1)/binning_factor-1))
                #print(int((i_index)/binning_factor)-1, int((j_index)/binning_factor)-1)
                binned_data[int((i_index)/binning_factor)-1, int((j_index)/binning_factor)-1] = binned_col

    return binned_data


def make_modes(data_file, n_modes, modes_out, image_shape=400, mode_shape=100):


    im = Image.open(data_file)
    im = im.convert('L')
    template = np.asarray(im, dtype=float)
    original = template/np.max(template)

    delta_deg = 0.72
    modes = np.zeros((50, 50, n_modes))
    n = 0
    heads = 2
    while n < n_modes - 1:
        head_size = int(mode_shape/heads)
        template = apply_binning(original, image_shape, head_size)
        super_head = np.zeros(((heads*2)*head_size, (heads*2)*head_size))
        size = head_size*heads
        half = mode_shape
        for i in range(heads*2):
            for j in range(heads*2):
                super_head[head_size*i:head_size*(i+1), head_size*j:head_size*(j+1)] = (-1)**i*template
        print(np.shape(super_head), head_size, size, half)
        im = Image.fromarray(super_head)
        for deg in np.arange(0, 360, delta_deg):
            rotated_im = im.rotate(deg)
            piece = np.asarray(rotated_im)[int(half-size/2):int(half+size/2), int(half-size/2):int(half+size/2)]
            mode_out = np.zeros((50, 50))
            mode_out[:np.shape(piece)[0], :np.shape(piece)[1]] = apply_binning(piece, 100, 50) 
            modes[:, :, n] = mode_out
            print(n)
            n += 1
            if n == n_modes - 1:
                break
        heads *= 2

    hdu = fits.HDUList(fits.PrimaryHDU(modes))
    hdu.writeto(modes_out)


files = glob.glob('*.jpg') + glob.glob('*.png') + glob.glob('*.jpeg')
files = ['stelter.jpg']
for infile in files:
    name = infile.split('.')[0]
    print(f'Working on the {name} basis.')
    modes_out = f'{name}_modes.fits'
    make_modes(infile, 2500, modes_out, mode_shape=100)
    print(f'Modes written to {modes_out}')
