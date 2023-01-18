from .read_families_from_file import get_families
from .aux_funs_for_runfile import get_files_path, load_jsons, check_num_parents, \
    convert_data_to_Als, remove_duplicate_children, insert_parents_data_to_their_haps, do_planB, write_binaries_to_file
from .validation import is_valid_family
from .dualHaplotype import DualHaplotype
from .no_parents import insert_data_to_single_allele_in_par_haps,  \
    associate_children_to_par_haps_while_parents_dont_exist
from .parents_exist import associate_children_to_par_haps_while_parents_exist
from .adding_children_data import add_child_data
from .haplotypes_to_Glstrign import create_gl_and_write_to_file


def run_GRAMM(input_path, out_files_path, alleles_names, races_dict, open_ambiguity_sim, is_serology):
    out_GLstr, out_binary = get_files_path(out_files_path)
    low2high, antigen2group, group2antigen = load_jsons()

    aux_tools = {'problematic_child': None, 'is_serology': is_serology, 'low2high': low2high,
                 'antigen2group': antigen2group, 'group2antigen': group2antigen, 'amb': {}, 'binary_dict': {},
                 'parent_has_empty_hap': {}}
    errors_in_families = {}

    families_dict = get_families(input_path, aux_tools)  # get families from input file

    # sorted keys of families_dict, for process the families in the order (dict is not ordered)
    sorted_keys = sorted(families_dict.keys(), key=lambda my_key: int(my_key))  # check that it sorted

    for idx_fam in sorted_keys:

        family = families_dict[idx_fam]  # get one family
        par_num = check_num_parents(family)
        aux_tools['problematic_child'] = None  # initialize for the new family

        # convert the data structure of each pair alleles_names from list to Als
        # 'Als' is a class which inherited from 'list', with adjusted attributes for alleles_names
        convert_data_to_Als(family)

        remove_duplicate_children(family, alleles_names)

        valid_family = is_valid_family(family, alleles_names, par_num, idx_fam, aux_tools, errors_in_families)
        if not valid_family:
            continue  # we do not want to analyze this family, so continue to the next family

        # create dual haplotype for father and mother (without data, yet)
        hapF, hapM = DualHaplotype(alleles_names), DualHaplotype(alleles_names)

        # ------------ insert data to parents haplotypes ------------
        if par_num > 0:
            insert_parents_data_to_their_haps(family, hapF, hapM, par_num)
        else:  # num_par == 0
            success_insert_single_allele, allele_inserted_to_parents = \
                insert_data_to_single_allele_in_par_haps(family, alleles_names, hapF, hapM, idx_fam, errors_in_families)
            if not success_insert_single_allele:
                continue
        # ------------ insert data to parents haplotypes ------------/

        # ---- go over children, associate with the haplotypes they inherited and add data to parents haplotypes ----
        children_keys = [k for k in family.keys() if k not in ['F', 'M']]
        # planB is used for insertion the children with partial associating with parents (1 or less parent).
        # after the loop over all the children, we will go over the children in planB and try to associate (and
        # add data) again, in hope that the parents now contain more data, that will help to associate those children
        planB = []
        error_in_adding_so_continue_to_next_family = False

        for child in children_keys:
            if par_num > 0:
                success_associate, associating = \
                    associate_children_to_par_haps_while_parents_exist(hapF, hapM, family[child], alleles_names)
            else:  # num_par == 0
                success_associate, associating = \
                    associate_children_to_par_haps_while_parents_dont_exist(hapF, hapM, family[child],
                                                                            allele_inserted_to_parents)

            if not success_associate:  # partial associating: with one or less parent
                planB.append(child)

            # after we know which haplotypes the _child_ inherited, we try to add data to them,
            # based on the parents and the _child_ data in these haplotypes
            success_adding = add_child_data(hapF, hapM, family[child], child, len(children_keys), associating, alleles_names, aux_tools,
                                            idx_fam, errors_in_families)
            if not success_adding:
                error_in_adding_so_continue_to_next_family = True
                break
        if error_in_adding_so_continue_to_next_family:
            continue

        for child in planB:  # explained above (in definition of 'planB')
            success_adding = do_planB(hapF, hapM, family[child], child, len(children_keys), alleles_names, aux_tools, idx_fam,
                                      errors_in_families)
            if not success_adding:
                error_in_adding_so_continue_to_next_family = True
                break
        if error_in_adding_so_continue_to_next_family:
            continue
        # ---- go over children, associate with the haplotypes they inherited and add data to parents haplotypes ----/

        # ---------- convert the data in the parents haplotypes to gl string and write to a file ----------
        create_gl_and_write_to_file(hapF, hapM, alleles_names, idx_fam, par_num, aux_tools, out_GLstr,
                                                        races_dict, open_ambiguity_sim, errors_in_families)
        # ---------- convert the data in the parents haplotypes to gl string and write to a file ----------/

    # load binary dict to a file
    write_binaries_to_file(out_binary, aux_tools)

    return errors_in_families, aux_tools




