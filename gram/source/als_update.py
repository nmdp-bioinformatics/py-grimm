from .general_aux_funs import equal_als


class Als(list):
    """
    This class illustrates a list of alleles_names.
    It's a different class because the alleles_names are required difference treatment:
    Two alleles_names considered as equal if they are not contradictory
    For example, 02 and 02:01 are equal. 03:01 and 03:02 are not.
    It causes to different implementations of some of the functions
    For example:
        "contain": checking if allele exists in alleles_names list (e.g. 02:01 exists in [03, 02])
        "sublist": checking if one list is a sublist of another
        and more
    In additional, there are more methods we want to add.
    Therefore, we build a class that inherited from "List", but adds and overrides some of the function
    """
    def __init__(self):
        super().__init__()

    def __contains__(self, value):
        """
        check if self contains value (override "in")
        :param value: value
        :return: true if contain, false otherwise
        """
        if value == "":  # empty value considers as contained in each Als list
            return True
        lst1 = [equal_als(ele, value) for ele in self]  # check equality with every item in Als
        return True if any(lst1) else False  # return True if equal to one at least

    def __eq__(self, other):
        """
        check equality (override "==")
        :param other: other Als list
        :return: true if equal, false otherwise
        """
        if self is None or other is None:
            return False
        if len(self) == len(other):  # a necessary condition for equality
            lst1 = [ele in other for ele in self]
            lst2 = [ele in self for ele in other]
            return True if all(lst1 + lst2) else False  # all item in self are in other, and vice versa
        return False

    def __ne__(self, other):
        """
        override "!="
        :param other: other Als
        :return: true if not equal, false else
        """
        return not (self == other)

    def __add__(self, other):
        """
        override "+". need to be override, cause otherwise created a list, not an Als
        :param other: other Als
        :return: new Als, with elements from self and other
        """
        new = Als()
        for ele in self:
            new.append(ele)
        for ele in other:
            new.append(ele)
        return new

    def empty_Als(self):  # (exists in "add_data_by_children" in lines 45, 47)
        """
        check if the Als is empty.
        (use it instead of the usual ways of checking empty list (len(obj) == 0, or: if not obj), for cases that
        we have ["", ""], and we want to consider it as empty allele)
        :return: true if empty, false else
        """
        for ele in self:
            if ele:
                return False
        return True

    def copy_a(self):
        """
        create a copy to an Als.
        use this function (and not use regular copy()) because we want that an Als() will be created, not a list
        :return: copied Als list
        """
        copy_l = Als()
        for ele in self:
            copy_l.append(ele)
        return copy_l

    def sub_lst(self, other):
        """
        check if self is a sub list of other
        :param other: other Als list
        :return: true if self is sub list of other
        """
        lst1 = [ele in other for ele in self]
        return True if all(lst1) else False  # return True if all the items in self are in other

    def index_a(self, value):
        """
        return the item index, if exists. if not exist: return -1
        :param value: value we want its index
        :return: the value index if exist. otherwise, -1
        """
        for i in range(len(self)):
            if equal_als(self[i], value):
                return i
        return -1

    def remove_a(self, value):
        """
        check if a value exists in als list, and if exist: remove from the list.
        pay attention that it does not return error if the value is not exist, unlike remove() in list
        :param value: value to remove
        """
        if self.index_a(value) != -1:
            self.remove(self[self.index_a(value)])
            # the complexity of the remove action is because of cases that the allele we want to remove and the
            # allele that exists in 'self' are equal in the meaning of alleles (like '02' and '02:01'), but not in the
            # regular meaning

    def remove_equal(self, value):
        """
        same to "remove_a" but ensure that values are identical.
        for example, in "remove_a", the values 01, 01:02 are equal, and here they are not.
        :param value: value to remove
        """
        idx_equal = -1
        for i in range(len(self)):
            if self[i] == value:
                idx_equal = i
        if idx_equal != -1:
            self.remove(self[idx_equal])

    def merge(self, other):
        """
        merge 'other' to 'self', with no repetitions.
        for example, self = [01:02, 04], other = [01, 03:01], so return [01:02, 04, 03:01]
        :param other: Als
        :return: merged Als
        """

        lst = self.copy_a()  # create a copy for not changing 'self'

        if not any(other):  # other = ["", ""], so we dont have anything to add
            return self

        inter, _ = self.intersection(other)
        if inter == 2:  # self = [01, 02, 03, 01], other = [01, 01]
            return self

        # 'other' has two identical values, so it's a special case
        if other[0] == other[1] and other[0]:  # second condition for check that the values in 'other' are not empty
            if other[0] in self:
                lst.append(other[0])  # self: [01, 02], other: [02, 02], so lst: [01, 02, 02]
            else:
                lst.extend(other)  # self: [01, 02], other: [03, 03], so lst: [01, 02, 03, 03]
            return lst

        for item in other:
            if item not in lst and item:  # second condition for not add empty value ("") to the merged lst
                lst.append(item)
        return lst

    def intersection(self, other):
        """
        find intersection between two Als. return the count of intersections, and the index
        (if there are some, so return the last).
        Note: we go over on 'other' and compare with 'self', so if other=[01,01], self=[01, 02], so intersections=2,
        but if other=[01,02], self=[01,01] so intersections=1
        :param other: Als
        :return: intersections count, index of the intersection
        """
        intersection_count = 0
        idx_intersection_in_other = None

        # if 'self' or 'other' (or they both) are empty (["", ""]), consider as no intersection
        if (not any(self)) or (not any(other)):
            return intersection_count, idx_intersection_in_other

        for allele in other:
            if allele != '' and allele in self:  # added: allele != ''
                intersection_count += 1
                idx_intersection_in_other = other.index_a(allele)  # if 2 intersections, idx will be the second
        return intersection_count, idx_intersection_in_other

    def count_a(self, value):
        """
        count how many times value appears in self
        :param value: value
        :return: counts times
        """
        count = 0
        for element in self:
            if equal_als(value, element):
                count += 1
        return count


