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
