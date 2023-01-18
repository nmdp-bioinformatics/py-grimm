# import matplotlib.pyplot as plt


def is_full_haplo(hap, locus):
    if all(loci in hap for loci in locus):
        return True
    return False


"""def Plot(values, titel, ylabel, xlabel, save_path):
        #r2 = np.array(prob)
        # a = np.hstack(r2)
        #plt.hist(r2, bins='auto')  # arguments are passed to np.histogram
        plt.plot(values)
        plt.title(titel)
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        #plt.xlim([0, 5])
        #plt.show()
        plt.savefig(save_path)
        plt.clf()"""


# cutt the dict to new_size
def cut_dict(dict_val, new_size=100000000):
    cut_dict = {}
    dict_val = sorted(dict_val.items(), key=lambda x: x[1], reverse=True)
    # dict_val.clear()
    while len(cut_dict) < new_size:
        cut_dict[dict_val[0][0]] = dict_val[0][1]
        del dict_val[0]
    """for key, value in cut_dict[key] = value:
        cut_dict[key] = value
        #del dict_val[0]
        if len(cut_dict) > new_size:
            break"""
    # dict_val = cut_dict
    return cut_dict
