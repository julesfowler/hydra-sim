Many Headed Hydra Modal Basis
-----------------------------

Check out the paper for more info: https://arxiv.org/pdf/2204.00575.pdf



If you want to run this code as is -- no frills -- on my headshot, you can do
that with : 

```
>>> python hydra.py 
``` 

If you want to try it on your own headshot, give the code an image of a typical
jpg/jpeg format cropped to 400x400 pixels:

```python
# import hyrda
from hydra import make_modes, test_basis

# Read in an image -- 400 x 400 
image_infile = 'path/to/you.jpg'
name = image_infile.split('.')[0]
print(f'Working on the {name} basis.')
modes_out = f'{name}_modes.fits'

# Make me a basis set with 2500 modes
make_modes(image_infile, 2500, modes_out)
print(f'Modes written to {modes_out}')
    
# Test it
turbulence = 'one_scene_turbulence_144x.fits'
approx, difference_image, percent_difference = test_basis(turbulence, modes_out)
print(f'The {name} hydra modal basis has a {percent_difference} from the original turbulence image.')
```
