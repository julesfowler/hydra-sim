# Given img 224x224 img scaled [-1,1]

mode_map = {}
mode_images = {}
mode = 0
full = 32
for i in [1, 2, 4, 8]:
    perms = list(itertools.product("01", repeat=i))
    for index in range(0, len(perms), i):
        print(index, index+int(i))
        cols = perms[index:index+int(i)]
        mode_map[mode] = np.array(cols)
        size = int(full/i)
        binlet = apply_binning(cropped_im, 224, size)
        mode_img = np.zeros((full,full))
        for index_i, col in enumerate(cols):
            for index_j, elem in enumerate(col):
                print(elem)
                if int(elem[0]) == 0:
                    mode_img[index_i*size:(index_i+1)*size, index_j*size:(index_j+1)*size] = binlet
                else:
                    mode_img[index_i*size:(index_i+1)*size, index_j*size:(index_j+1)*size] = binlet.transpose()
        mode_images[mode] = mode_img.reshape(full,full)
        plt.clf()
        plt.imshow(mode_img.reshape(full,full))
        plt.savefig(f'mode_{mode}.png')
        mode += 1


