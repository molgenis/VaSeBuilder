class SuperContext:
    def __init__(self, first_varcon, second_varcon):
        """Constructor that creates the super context from two variant contexts.

        The initial super context data such as chromomsome name, origin, start and end position will be copied from the
        first variant context, which is then also linked to the super context. The second variant context will be
        linked to the super context as afterwards and may overwrite the start and/or end position of the merged context.

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
        """Returns the chromosome name of the super context.

        Returns
        -------
        self.super_chrom : str
            Chromosome name
        """
        return self.super_chrom

    def get_origin(self):
        """Returns the origin position of the super context.

        The origin position is based on the origin position of the first variant context, that in essence created the
        merged context.

        Returns
        -------
        self.super_origin : int
            Origin point of the first variant context
        """
        return self.super_origin

    def get_start(self):
        """Returns the leftmost genomic position of the super context.

        Returns
        -------
        self.super_start : int
            Super context leftmost position
        """
        return self.super_start

    def get_end(self):
        """Returns the rightmost genomic position of the super context.

        Returns
        -------
        self.super_end : int
            Super context rightmost position
        """
        return self.super_end

    def get_variant_contexts(self):
        """Returns all variant contexts associated with the super context.

        Returns
        -------
        self.variant_contexts: list of VariantContext
            Variant contexts associated with the super context
        """
        return self.variant_contexts

    def get_num_of_variant_contexts(self):
        """Returns the number of variant contexts associated with the current super context.

        Returns
        -------
        int
            The number of associated variant contexts
        """
        return len(self.variant_contexts)

    def get_context(self):
        """Returns a context array of the super context.

        Returns
        -------
        list of str and int
            Returns a context array of the super context
        """
        return [self.super_chrom, self.super_origin, self.super_start, self.super_end]

    def get_genomic_region(self):
        """Returns the super context as a genomic region representation string.

        The genomic region string consists of chrom:start-end (i.e.: '21:100-1000')

        Returns
        -------
        str
            Genomic region representation (i.e.: 21:100-1000)
        """
        return f"{self.super_chrom}:{self.super_start}-{self.super_end}"

    def add_variant_context(self, varcon_toadd):
        """Adds a variant context to the current super context.

        Prior to adding, it is checked whether the variant context is on the same chromosome as the super context.
        Furthermore, the new super context start and end positions are also set if the start and end of the variant
        context are smaller and larger than the super context start and end respectively.

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
        """Adds a list of variant contexts

        Parameters
        ----------
        varcons_toadd : list of VariantContexts
        """
        for varcon in varcons_toadd:
            self.add_variant_context(varcon)

    def determine_super_start(self, new_start):
        """Determines the new super context start position.

        A new start position is only set if the provided start position is smaller than the current super context start
        position.

        Parameters
        ----------
        new_start : int
            New start position to check
        """
        if new_start < self.super_start:
            self.super_start = new_start

    def determine_super_end(self, new_end):
        """Determine the new super context end position.

        A new end position is only set if the provided end position is larger than the current super context end
        position.

        Parameters
        ----------
        new_end : int
            New end position to check
        """
        if new_end > self.super_end:
            self.super_end = new_end

    def remove_variant_context(self, varcon_toremove):
        """Removes a specified variant context from the super context.

        Parameters
        ----------
        varcon_toremove: VariantContext
            Variant context to remove
        """
        print("")

    def split_super_context(self, varcon_touse, keep_left=True):
        """Splits the super context into two via a provided variant context.

        The provided variant context is used as a breaking point to split the super context into two smaller super
        contexts. If the parameter keep_varcon is set to True, the variant context used as a breaking point is kept in
        the current super context. If set to false, than the variant context used as a breaking point is placed in the
        other super context.

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

        # Check whether the starting variant context should be added to the left or right contexts
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
            split_context = self.construct_super_context(right_contexts[0], right_contexts[1], right_contexts[2:])
        else:
            right_contexts[0].remove_from_supercontext()
            self.remove_variant_context(right_contexts[0])

        return split_context

    def get_left_contexts(self, start_varcon):
        """Returns all variant contexts left from the specified variant context.

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
        """Returns all variant context right from the specified variant context.

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

    def construct_super_context(self, first_varcon, second_varcon, remaining_varcons):
        """Constructs and returns a super context.

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

    def fix_gaps(self):
        """Checks if there is a gap in the super context.

        If the super context has a gap, that can be introduced when a variant context is removed, the super context is
        split using a similar method to split_super_context().

        Returns
        -------
        SuperContext or None
        """

    def has_gaps(self):
        """Checks and returns whether the super context has gaps.

        Returns
        -------
        bool
            True if the super context has gaps, False if not
        """
        left_positions = []
        varcon_per_leftpos = {}

        # First gather all variant contexts by leftmost genomic position
        for varcon in self.variant_contexts:
            vcleftpos = varcon.get_variant_context_start()
            left_positions.append(vcleftpos)
            if vcleftpos not in varcon_per_leftpos:
                varcon_per_leftpos[vcleftpos] = []
            varcon_per_leftpos[vcleftpos].append(varcon)
        left_positions = list(set(left_positions))

        # Start checking for gaps by
        for x in range(len(left_positions)):
            # Get the left context
            left_varcon = varcon_per_leftpos[left_positions[x]]
            if len(left_varcon) > 1:
                left_varcon = self.determines_largest_gap_context(left_varcon)

            # Get the right context of the current list position
            right_varcon = varcon_per_leftpos[left_positions[x+1]]

    def determines_largest_gap_context(self, varconlist):
        return "aap"
