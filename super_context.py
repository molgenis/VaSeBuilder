"""WIP multi-context."""



class SuperContext:
    """WIP multi-context."""

    def __init__(self, first_varcon, second_varcon):
        """Create the super context from two variant contexts.

        The initial super context data such as chromosome name, origin, start
        and end position will be copied from the first variant context, which
        is then also linked to the super context. The second variant context
        will be linked to the super context as afterwards and may overwrite
        the start and/or end position of the merged context.

        Parameters
        ----------
        first_varcon : VariantContext
            First variant context to add
        second_varcon: VariantContext
            Second variant context to add
        """
        self.variant_contexts = []
        self.super_chrom = first_varcon.get_variant_context_chrom()
        self.super_origin = first_varcon.get_variant_context_origin()
        self.super_start = first_varcon.get_variant_context_start()
        self.super_end = first_varcon.get_variant_context_end()
        self.add_variant_context(second_varcon)

    def get_chrom(self):
        """Return the chromosome name of the super context.

        Returns
        -------
        self.super_chrom : str
            Chromosome name
        """
        return self.super_chrom

    def get_origin(self):
        """Return the origin position of the super context.

        The origin position is based on the origin position of the first
        variant context, that in essence created the merged context.

        Returns
        -------
        self.super_origin : int
            Origin point of the first variant context
        """
        return self.super_origin

    def get_start(self):
        """Return the leftmost genomic position of the super context.

        Returns
        -------
        self.super_start : int
            Super context leftmost position
        """
        return self.super_start

    def get_end(self):
        """Return the rightmost genomic position of the super context.

        Returns
        -------
        self.super_end : int
            Super context rightmost position
        """
        return self.super_end

    def get_length(self):
        """Determine and return the length of the super context.

        Returns
        -------
        int
            Length of the super context
        """
        return self.super_end - self.super_start

    def get_variant_contexts(self):
        """Return all variant contexts associated with the super context.

        Returns
        -------
        self.variant_contexts: list of VariantContext
            Variant contexts associated with the super context
        """
        return self.variant_contexts

    def get_num_of_variant_contexts(self):
        """Return the number of variant contexts associated with the current super context.

        Returns
        -------
        int
            The number of associated variant contexts
        """
        return len(self.variant_contexts)

    def get_context(self):
        """Return a context array of the super context.

        Returns
        -------
        list of str and int
            Returns a context array of the super context
        """
        return [self.super_chrom, self.super_origin, self.super_start, self.super_end]

    def get_genomic_region(self):
        """Return the super context as a genomic region representation string.

        The genomic region string consists of chrom:start-end (i.e.: '21:100-1000')

        Returns
        -------
        str
            Genomic region representation (i.e.: 21:100-1000)
        """
        return f"{self.super_chrom}:{self.super_start}-{self.super_end}"

    def add_variant_context(self, varcon_toadd):
        """Add a variant context to the current super context.

        Prior to adding, it is checked whether the variant context is on the
        same chromosome as the super context. Furthermore, the new super
        context start and end positions are also set if the start and end of
        the variant context are smaller and larger than the super context
        start and end respectively.

        Parameters
        ----------
        varcon_toadd : VariantContext
            Variant context to add to the super context
        """
        if varcon_toadd.get_variant_context_chrom() == self.super_chrom:
            self.variant_contexts.append(varcon_toadd)
            self.determine_super_start(varcon_toadd.get_variant_context_start())
            self.determine_super_end(varcon_toadd.get_variant_context_end())

    def add_variant_contexts(self, varcons_toadd):
        """Add a list of variant contexts.

        Parameters
        ----------
        varcons_toadd : list of VariantContexts
            Variant contexts to add to the current super context
        """
        for varcon in varcons_toadd:
            self.add_variant_context(varcon)

    def determine_super_start(self, new_start):
        """Determine the new super context start position.

        A new start position is only set if the provided start position is
        smaller than the current super context start position.

        Parameters
        ----------
        new_start : int
            New start position to check
        """
        if new_start < self.super_start:
            self.super_start = new_start

    def determine_super_end(self, new_end):
        """Determine the new super context end position.

        A new end position is only set if the provided end position is larger
        than the current super context end position.

        Parameters
        ----------
        new_end : int
            New end position to check
        """
        if new_end > self.super_end:
            self.super_end = new_end

    # TODO:
    def remove_variant_context(self, varcon_toremove):
        """Remove a specified variant context from the super context.

        Parameters
        ----------
        varcon_toremove: VariantContext
            Variant context to remove from the super context
        """
        print("thing to remove")

    def split_super_context(self, varcon_touse, keep_left=True):
        """Split the super context into two via a provided variant context.

        The provided variant context is used as a breaking point to split the
        super context into two smaller super contexts. If the parameter
        keep_varcon is set to True, the variant context used as a breaking
        point is kept in the current super context. If set to false, than the
        variant context used as a breaking point is placed in the other
        super context.

        Parameters
        ----------
        varcon_touse : VariantContext
            Variant context to use for splitting
        keep_left : bool
            Whether to keep the variant context in the current super context

        Returns
        -------
        split_context : SuperContext or None
            A super context that was split from the current super context
        """
        left_contexts = self.get_left_contexts(varcon_touse)
        right_contexts = self.get_right_contexts(varcon_touse)
        split_context = None

        # Check whether the starting variant context should be added to the
        # left or right contexts
        if keep_left:
            left_contexts.append(varcon_touse)
        else:
            right_contexts.insert(0, varcon_touse)

        # Determine what to do with the left contexts
        if len(left_contexts) >= 2:
            self.variant_contexts = left_contexts
        else:
            left_contexts[0].remove_from_supercontext()
            self.remove_variant_context(left_contexts[0])
            self.variant_contexts = []

        # Determine what to do with the right contexts
        if len(right_contexts) >= 2:
            split_context = self.construct_super_context(right_contexts[0],
                                                         right_contexts[1],
                                                         right_contexts[2:])
        else:
            right_contexts[0].remove_from_supercontext()
            self.remove_variant_context(right_contexts[0])

        return split_context

    def get_left_contexts(self, start_varcon):
        """Return all variant contexts left from the specified variant context.

        Parameters
        ----------
        start_varcon : VariantContext
            Variant context as starting point

        Returns
        -------
        left_contexts : list of VariantContexts
            Variant contexts to the left
        """
        svend = start_varcon.get_variant_context_end()
        left_contexts = []
        for varcon in self.variant_contexts:
            if varcon.get_variant_context_end() <= svend:
                left_contexts.append(varcon)
        return left_contexts

    def get_right_contexts(self, start_varcon):
        """Return all variant context right from the specified variant context.

        Parameters
        ----------
        start_varcon : VariantContext
            Variant context as starting point

        Returns
        -------
        right_contexts : list of VariantContext
            Variant contexts to the right
        """
        svend = start_varcon.get_variant_context_end()
        right_contexts = []
        for varcon in self.variant_contexts:
            if varcon.get_variant_context_start() <= svend:
                right_contexts.append(varcon)
        return right_contexts

    @staticmethod
    def construct_super_context(first_varcon, second_varcon, remaining_varcons):
        """Construct and returns a super context.

        Parameters
        ----------
        first_varcon : VariantContext
            First variant context to create super context
        second_varcon : VariantContext
            Second variant context
        remaining_varcons : list of VariantContext or None
            Remaining variant contexts to add

        Returns
        -------
        super_context : SuperContext
            Newly constructed super context
        """
        super_context = SuperContext(first_varcon, second_varcon)
        if remaining_varcons is not None:
            if len(remaining_varcons) > 0:
                super_context.add_variant_contexts(remaining_varcons)
        return super_context

    def has_gaps(self):
        """Check and return whether the current merged context has a gap.

        Variant contexts associated with the merged context are first sorted
        on leftmost genomic position. Then, overlap between the current
        variant context and its first neighbor on the right. The method
        stops after the first gap has been found as that indicates the gap
        or gaps need to be fixed.

        Returns
        -------
        bool
            True if a gap has been found, False if not
        """
        # Gather the left most positions and the variant contexts per leftmost position
        vcleft_varcons = self.get_varcons_by_leftpos()
        vcleftpos = list(set(vcleft_varcons.keys()))

        # Start determining whether there are any gaps
        curvarcon = vcleft_varcons[vcleftpos[0]]
        curvarcon = self.determine_largest_variantcontext(curvarcon)

        for vclp in vcleftpos[1:]:
            varcon = vcleft_varcons[vclp]
            varcon = self.determine_largest_variantcontext(varcon)

            # Check for overlap between curvarcon and the selected varcon.
            if self.varcons_overlap(curvarcon, varcon):
                curvarcon = self.determine_rightmost_context(curvarcon, varcon)
            else:
                return True
        return False

# TODO: =======================================================================
#     def determines_largest_gap_context(self, varconlist):
#         return "test"
# =============================================================================

    @staticmethod
    def select_varcons_for_super_context(leftpositions, varcons_by_leftpos):
        """Return variant contexts to create a new super context with.

        Parameters
        ----------
        leftpositions : list of int
            Leftmost genomic positions from variant contexts
        varcons_by_leftpos : dict
            Variant contexts to construct super context with

        Returns
        -------
        varcons_for_sc : list of VariantContext
            Variant contexts to put in a SuperContext
        """
        varcons_for_sc = []
        for leftpos in leftpositions:
            varcons_for_sc.extend(varcons_by_leftpos[leftpos])
        return varcons_for_sc

    @staticmethod
    def determine_rightmost_context(first_varcon, second_varcon):
        """Determine and return the variant context with the rightmost genomic position.

        Parameters
        ----------
        first_varcon : VariantContext
            First variant context to use
        second_varcon : VariantContext
            Second variant context to use

        Returns
        -------
        VariantContext
            Variant context that has the rightmost genomic position
        """
        if second_varcon.get_variant_context_end() > first_varcon.get_variant_context_start():
            return second_varcon
        return first_varcon

    @staticmethod
    def varcons_overlap(first_varcon, second_varcon):
        """Determine whether two provided variant contexts overlap.

        Parameters
        ----------
        first_varcon : VariantContext
            First variant context
        second_varcon : VariantContext
            Second variant context

        Returns
        -------
        bool
            True if variant contexts overlap, False if not
        """
        start1 = first_varcon.get_variant_context_start()
        stop1 = first_varcon.get_variant_context_end()
        start2 = second_varcon.get_variant_context_start()
        stop2 = second_varcon.get_variant_context_end()
        if start2 <= stop1 and start1 <= stop2:
            return True
        return False

    def determine_largest_variantcontext(self, variant_contexts):
        """Determine and return the largest variant context of a list of variant contexts.

        The largest context is determined based on the leftmost and rightmost
        genomic position. If there is only one variant context, that variant
        context is returned.

        Parameters
        ----------
        variant_contexts : list of VariantContext
            Variant contexts to obtain the largest context from

        Returns
        -------
        largest_context : VariantContext or None
            The largest variant context, None if the list is empty
        """
        largest_context = None
        if len(variant_contexts) == 1:
            largest_context = variant_contexts[0]
        elif len(variant_contexts) >= 2:
            largest_context = None
            for varcon in variant_contexts:
                if largest_context is not None:
                    largest_context = self.get_largest_context(largest_context, varcon)
                else:
                    largest_context = varcon
        return largest_context

    @staticmethod
    def get_largest_context(first_context, second_context):
        """Determine and return the largest of two variant contexts.

        Parameters
        ----------
        first_context : VariantContext
            First variant context to check
        second_context : VariantContext
            Other variant context to check

        Returns
        -------
        VariantContext
            Largest variant context of two variant contexts
        """
        if (second_context.get_variant_context_length() >=
                first_context.get_variant_context_length()):
            return second_context
        return first_context

    def get_varcons_by_leftpos(self):
        """Return variant contexts in the super context by leftmost genomic position.

        Returns
        -------
        varcons_by_leftpos : dict
            Varcons per leftmost genomic position
        """
        varcons_by_leftpos = {}
        for varcon in self.variant_contexts:
            varcon_leftpos = varcon.get_variant_context_start()
            if varcon_leftpos not in varcons_by_leftpos:
                varcons_by_leftpos[varcon_leftpos] = []
            varcons_by_leftpos[varcon_leftpos].append(varcon)
        return varcons_by_leftpos

    def fix_gaps(self):
        """Fix gaps in the current super context.

        Each gap in the current super context is fixed by splitting the super
        context in two smaller super contexts.

        Returns
        -------
        fixed_super_contexts : list of SuperContext
            List of SuperContexts if one or more gaps were found,
            empty list if no gaps were found
        """
        fixed_super_contexts = []
        leftindex = 0

        # Gather the left most positions and the variant contexts per leftmost position
        vcleft_varcons = self.get_varcons_by_leftpos()
        vcleftpos = list(set(vcleft_varcons.keys()))

        # Obtain the next variant context
        curvarcon = vcleft_varcons[vcleftpos[0]]
        curvarcon = self.determine_largest_variantcontext(curvarcon)

        # Start determining whether there are any gaps
        for vclp in vcleftpos[1:]:
            varcon = vcleft_varcons[vclp]
            varcon = self.determine_largest_variantcontext(varcon)

            # Check for overlap between curvarcon and the selected varcon
            if self.varcons_overlap(curvarcon, varcon):
                curvarcon = self.determine_rightmost_context(curvarcon, varcon)
            else:
                curindex = vcleftpos.index(vclp)
                varcons_for_super_context = self.select_varcons_for_super_context(
                    vcleftpos[leftindex:curindex],
                    vcleft_varcons
                    )

                # We found a gap and left context is alone
                if len(varcons_for_super_context) == 1:
                    # Remove variant context from the SuperContext it belonged to
                    varcons_for_super_context[0].remove_from_supercontext()
                    self.remove_variant_context(varcons_for_super_context[0])

                # We found a gap and need to create a new smaller super context
                elif len(varcons_for_super_context) >= 2:
                    fixed_super_contexts.append(SuperContext(varcons_for_super_context[0],
                                                             varcons_for_super_context[1]))

                # Update the required data
                curvarcon = varcon
                leftindex = vcleftpos.index(vclp)
        return fixed_super_contexts
