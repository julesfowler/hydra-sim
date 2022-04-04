## -- IMPORTS 
from astropy.io import fits
import glob
import math
import matplotlib.pyplot as plt
import numpy as np

from PIL import Image

## -- MAIN FUNCTIONS

def apply_binning(data, resolution, modes):
    """ Binning algorithm based on Maaike van Kooten's AO sumer school
    notebook demo. Bins data down to a given size. *NOTE* this function 
    trusts you to do what you've told it. It will not do any smart type 
    checking it will throw you an error or build a bad image if you don't
    give it the sizes you promised. 

    Parameters
    ----------
    data : np.array
        Input image to bin down, should be resolution x resolution. 
    resolution : int
        Original resolution / size of image. (Assumes square input image.)
    modes : int
        Final resolution / size of image. (Assumes square output image.)

    Returns 
    -------
    binned_data : np.array
        Output binned image binned modes x modes.
    """
    
    # Check to see if there's any binning to be done at all.
    if resolution==modes:
        binned_data = data

    else:
        
        # Figure out the super pixels for the new image
        binned_data = np.zeros((modes, modes))
        binning_factor = int(resolution/modes)
        
        # Run through and build super pixels
        for i_index in np.arange(binning_factor, resolution+1, binning_factor):
            for j_index in np.arange(binning_factor, resolution+1, binning_factor):
                binned_col = np.mean(
                             data[i_index-binning_factor:i_index,
                                  j_index-binning_factor:j_index]
                             )
                # Update each new pixel
                binned_data[int((i_index)/binning_factor)-1, int((j_index)/binning_factor)-1] = binned_col
    
    # Return new image
    return binned_data


def image_decomposition(original_image, A, A_i, thresh=1e-5):
    """ Function to decompose an image into a basis. Originally written for my
    real research, in julesfowler/powerseal. Solves Ax = b for x, given b is the
    image and A is the modal basis. 

    Parameters
    ----------
    original_image : np.array
        Image to decompose.
    A : np.array
        2D basis set. 
    A_i : np.array
        Inverse (or more likely pseudo-inverse) of A.
    thresh : float, optional
        Whatever threshold to trigger the warning. 
    
    Returns
    -------
    x : np.array
        Coefficients that make up x.
    approx :  np.array
        Approximate image built by the decomposition. 
    """
    # Flatten the image into a vector
    b = original_image.flatten()
    
    # Caclulate x with x = A_inverse*b
    x = np.matmul(A_i, b)
    
    # Calculate the approximation of the image with this decomposition
    approx = np.matmul(A, x)
    
    # Throw a warning if the image decomposition isn't very good. 
    # For a real basis set (i.e., Fourier/etc) this should be on the order of
    # machine precision. 
    if np.mean(np.abs(b - approx)) > thresh:
        print("Approximation isn't up to standard.")
        print(f'{np.mean(np.abs(b - approx))} > {thresh}')
    
    # Return the solution, and the approximate image
    return x, approx


def make_modes(data_file, n_modes, modes_out, image_shape=400, mode_shape=100, final_bin=50, smallest_head=2):
    """ Builds a modal basis given an input image file. *NOTE* this function
    trusts you to feed it a square image of the size you've promised. (In
    another life I would build more intelligent functions but this is an April
    Fools day package after all...) 

    Parameters
    ----------
    data_file : str
        Path to the image file. Needs a standard jpg/png/jpeg/pdf image type.
    n_modes : int
        The number of modes (along one dimension) to build.     
    modes_out : str
        Filename to write out mode data to. 
    image_shape : int, optional
        One dimension of the image, in pixels. Requires a square image. Defaults
        to 400.
    mode_shape : int
        One dimension of the size of the mode to build out, in pixels. (Note
        that if you apply a final bin, this will actually be the shape of the
        mode out.) Defaults to 100.
    final_bin : bool or int, optional
        Whether or not to apply a final binning to the modes, to make them
        smaller. Use this option to get smaller modes out of starting images
        that don't divide nicely. 
    smallest_head : int
        The smallest your recursive head can be. This really shouldn't be
        smaller than 2 pixels, you can make it bigger to preserve more
        structure. 
        
    Returns
    -------
    modes : np.array
        An (mode_shape, mode_shape, n_modes) sized array containing the modal
        basis set.
    """
    
    # Read in image, convert to black and white, cast as array, scale 0-1
    im = Image.open(data_file)
    im = im.convert('L')
    template = np.asarray(im, dtype=float)
    original = template/np.max(template)
    
    # Find out how many heads will be in the highest frequency image
    # Given the final mode cannot have heads smaller than smallest_head
    final_mode_shape = final_bin if final_bin else mode_shape
    high_frequency_heads = np.floor(math.log(final_mode_shape, smallest_head))

    # This sets the rotation frequency to build as many modes as you ask for
    delta_deg = 360/(n_modes/high_frequency_heads)
    
    # Start off building modes with 2 heads
    modes = np.zeros((final_mode_shape, final_mode_shape, n_modes))
    n = 0
    heads = 2
    
    # Run through until we've built enough modes 
    # Since we're using floor we may not complete the sequence before we have enough
    while n < n_modes - 1:
        
        # Figure out the head size for this iteration and bin down the head
        head_size = int(mode_shape/heads)
        template = apply_binning(original, image_shape, head_size)
        
        # Build a mode twice as big so we fill the whole area as we rotate
        super_head = np.zeros(((heads*2)*head_size, (heads*2)*head_size))
        size = head_size*heads
        half = mode_shape

        # Build up rows of head / inverse head
        for i in range(heads*2):
            for j in range(heads*2):
                super_head[head_size*i:head_size*(i+1), head_size*j:head_size*(j+1)] = (-1)**i*template
        
        # Convert back to Image type for ease of rotation
        im = Image.fromarray(super_head)

        # Fill out the frequency through each rotation incriment 
        for deg in np.arange(0, 360, delta_deg):
            rotated_im = im.rotate(deg)

            # Cut the mode down to its actual size
            piece = np.asarray(rotated_im)[int(half-size/2):int(half+size/2), int(half-size/2):int(half+size/2)]
            
            # Pad edges with zeros if there are edge effects
            mode_out = np.zeros((final_mode_shape, final_mode_shape))
            # Apply final bin if specified 
            mode_out[:np.shape(piece)[0], :np.shape(piece)[1]] = apply_binning(piece, 100, final_mode_shape) 
            modes[:, :, n] = mode_out
            n += 1
            
            # See if we've built enough modes
            if n == n_modes - 1:
                break
        heads *= 2
    
    # Remove any nans 
    modes[np.isnan(modes)] = 0.0

    # Write out the modes 
    hdu = fits.HDUList(fits.PrimaryHDU(modes))
    hdu.writeto(modes_out, overwrite=True)
    
    # Return them
    return modes


def test_basis(turbulence_image, mode_data, image_size=100, rcond=1e-6):
    """ Function to decompose an example turbulence image as a test of the basis
    set. 

    Parameters
    ----------
    turbulence_image : str
        Path to a turbulence image. 
    
    mode_data : str
        Path to a fits file with modes. 

    Returns 
    -------
    approx : np.array
        Appromxiation of the turbulence image. 
    difference_image : np.array
        Difference image between the approximation and the error. 
    percent_difference : float
        Percent difference of the approximation from the original.
    """
    
    # Pull out plotting name
    name = mode_data.split('_')[0]
    
    # Pull data from fits files
    modes = fits.getdata(mode_data)
    turbulence = fits.getdata(turbulence_image)
    mode_size = np.shape(modes)[0]
    binned_turbulence = apply_binning(turbulence[:image_size, :image_size], image_size, mode_size)
    A = modes.reshape(mode_size**2, mode_size**2)
    pinv = np.linalg.pinv(A, rcond=1e-6)
    x, approx = image_decomposition(binned_turbulence, A, pinv, thresh=1e-5)

    difference_image = approx.reshape(mode_size, mode_size) - binned_turbulence 
    percent_difference = 100*(np.median(difference_image)/np.median(binned_turbulence))
    
    return approx, difference_image, percent_difference 


if __name__ == "__main__":
    
    # Read me in 
    image_infile = 'fowler.jpeg'
    name = image_infile.split('.')[0]
    print(f'Working on the {name} basis.')
    modes_out = f'{name}_modes.fits'

    # Make me a basis set with 2500 modes
    make_modes(image_infile, 2500, modes_out)
    print(f'Modes written to {modes_out}')
    
    # Test me 
    turbulence = 'one_scene_turbulence_144x.fits'
    approx, difference_image, percent_difference = test_basis(turbulence, modes_out)
    print(f'The {name} hydra modal basis has a {percent_difference} from the original turbulence image.')
    
