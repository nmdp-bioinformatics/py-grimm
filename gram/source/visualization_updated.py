from matplotlib import pyplot as plt
import matplotlib.patches as patches


def visualize(resultsFile, parent_has_empty_hap, output_path):
    SOL_PATH = resultsFile
    # change to solution number (some families have many solutions)
    SOL_INDEX = 0

    SOL_FILE = open(SOL_PATH, 'r+')
    SOL_LINES = SOL_FILE.readlines()
    cur_FAMCODE = -1
    threshold_fig = 0
    results_exist = False
    for idx, sol_line in enumerate(SOL_LINES):
        if threshold_fig == 5:
            break  # only 5 figures
        if cur_FAMCODE == sol_line.split(',')[0]:
            continue  # visualization to the first result of each family only
        cur_FAMCODE = sol_line.split(',')[0]
        fam_solutions = []
        if sol_line.split(',')[0] != "FAMCODE":  # skip header line
            results_exist = True
            fam_solutions.append(sol_line.split(',')[1:])
            FAM_CODE = sol_line.split(',')[0]
            CHILDS_INDEX = 4
            childs_num = fam_solutions[SOL_INDEX][CHILDS_INDEX].count('=')

            # initialize graph components sizes
            mid = 0.5
            parent_offset = 0.05
            parent_width = 0.4
            height = 0.15
            parent_top = 0.8

            # initialize graph size
            plt.figure(idx, figsize=(17, 7))
            currentAxis = plt.gca()
            currentAxis.set_axis_off()

            # create rectangles for each parent and add haplotype text
            colors_f = ['#3D60A8', '#7EA0C5']
            colors_m = ['#A33F62', '#CEACB0']

            # if there is unknown haplotype, its color is white (#FFFFFF)
            if FAM_CODE in parent_has_empty_hap and parent_has_empty_hap[FAM_CODE] == 'F':
                colors_f[1] = '#FFFFFF'
            elif FAM_CODE in parent_has_empty_hap and parent_has_empty_hap[FAM_CODE] == 'M':
                colors_m[1] = '#FFFFFF'

            for hap in range(2):
                father_coords = (mid-parent_offset-parent_width,parent_top-height-hap*height)
                mother_coords = (mid+parent_offset,parent_top-height-hap*height)
                father_text_coords = ((father_coords)[0]+parent_width/2,father_coords[1]+height/2)
                mother_text_coords = (mother_coords[0]+parent_width/2,mother_coords[1]+height/2)
                if hap == 0:
                    father_caption = (father_text_coords[0],father_text_coords[1]+height)
                    mother_caption = (mother_text_coords[0],mother_text_coords[1]+height)
                    currentAxis.annotate('Father',xy=father_caption,color='black', weight='bold', fontsize=10, ha='center',va='center')
                    currentAxis.annotate('Mother',xy=mother_caption,color='black', weight='bold', fontsize=10, ha='center',va='center')
                father_hap_rect = patches.Rectangle(father_coords,parent_width,height,linewidth=1,edgecolor='black',facecolor=colors_f[hap])
                currentAxis.annotate('F'+str(hap+1)+': '+fam_solutions[SOL_INDEX][hap].strip('"'),xy=father_text_coords,color='black', weight='bold', fontsize=10, ha='center',va='center')
                mother_hap_rect = patches.Rectangle(mother_coords,parent_width,height,linewidth=1,edgecolor='black',facecolor=colors_m[hap])
                currentAxis.annotate('M'+str(hap+1)+': '+fam_solutions[SOL_INDEX][2+hap].strip('"'),xy=mother_text_coords,color='black', weight='bold', fontsize=10, ha='center',va='center')
                #ax.text(father_coords[0],father_coords[1],'Hello',fontsize=24)
                currentAxis.add_patch(father_hap_rect)
                currentAxis.add_patch(mother_hap_rect)

            # create rectangle for each _child_ and add haplotype text
            childs_right = 0.1
            childs_space = 0.02
            childs_top = 0.35
            child_width = 0.15

            if childs_num > 5:
                childs_right = 0.1
                childs_space = 0.02
                childs_top = 0.35
                child_width = 0.09

            dict_colors = {'F1': colors_f[0], 'F2': colors_f[1], 'M1': colors_m[0], 'M2': colors_m[1]}
            for child_num in range(childs_num):
                cur_child = fam_solutions[SOL_INDEX][CHILDS_INDEX].split(';')[child_num].lstrip('"')
                for hap in range(2):
                    child_hap_fm = cur_child.split("=")[1].split("~")[hap].rstrip()
                    try:
                        child_color = dict_colors[child_hap_fm]
                    except:
                        child_color = 'none'
                    # complex computation of x,y for current _child_
                    child_coords = (mid+(childs_space+child_width)*(child_num-childs_num/2+0.5)-child_width/2,childs_top-height-hap*height)
                    child_text_coords = (child_coords[0]+child_width/2,child_coords[1]+height/2)
                    if hap == 0:
                        child_caption = (child_text_coords[0],child_text_coords[1]+height)
                        currentAxis.annotate('Child '+cur_child.split('=')[0].lstrip('C'),xy=child_caption,color='black', weight='bold', fontsize=10, ha='center',va='center')
                    child_hap_rect = patches.Rectangle(child_coords,child_width,height,linewidth=1,edgecolor='black', facecolor=child_color)
                    hap_text = cur_child.split('=')[1].split('~')[hap].replace('+'," and ").rstrip('\n').rstrip('"')
                    currentAxis.add_patch(child_hap_rect)
                    currentAxis.annotate(hap_text,xy=child_text_coords,color='black', weight='bold', fontsize=10, ha='center',va='center')
            # output mark
            # ax.annotate("Output",xy=(0,0.5),color='black', weight='bold', fontsize=24, ha='center',va='center')

            threshold_fig += 1
            name_png = output_path + '/' + FAM_CODE + '.png'
            name_pdf = output_path + '/' + FAM_CODE + '.pdf'
            plt.savefig(name_png, dpi=300)
            plt.savefig(name_pdf)
            plt.close('all')
            # plt.show()

    # no results, so create png with the massage: 'No data to show'
    if not results_exist:
        # initialize graph size
        plt.figure(idx, figsize=(17, 7))
        currentAxis = plt.gca()
        currentAxis.set_axis_off()
        currentAxis.annotate('No data to show', xy=(0.5, 0.5), color='black', weight='bold', fontsize=20, ha='center',
                             va='center')
        name = output_path + '/1.png'
        plt.savefig(name, dpi=300)
        plt.close('all')
