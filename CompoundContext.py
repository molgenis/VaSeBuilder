class CompoundContext:
    def __init__(self, first_varcon, second_varcon):
        self.variant_contexts = [first_varcon, second_varcon]
        self.priority_level = 64
        self.priority_label = ""

    def get_variant_contexts(self):
        """Returns the variant contexts forming a compound context.

        Returns
        -------
        self.variant_contexts : list of VariantContext
            Variant contexts forming a
        """
        return self.variant_contexts

    def get_num_of_variant_contexts(self):
        """Returns the number of variant contexts forming a compound context.

        Returns
        -------
        int
            The number of variant contexts forming a compound context
        """
        return len(self.variant_contexts)

    def add_variant_context(self, varcon_toadd):
        """Adds a provided VariantContext to the compound context.

        Parameters
        ----------
        varcon_toadd : VariantContext
            self.variant_contexts.append(varcon_toadd)
        """
        if varcon_toadd is not None:
            self.variant_contexts.append(varcon_toadd)

    def add_variant_contexts(self, varcons_to_add):
        """Adds a list of variant contexts to a compound context.

        Parameters
        ----------
        varcons_to_add : list of VariantContext
            List of variant contexts to add to the compound
        """
        if varcons_to_add is not None:
            if len(varcons_to_add) >= 1:
                self.variant_contexts.extend(varcons_to_add)

    def determine_priority_level(self):
        print("aap")

    def remove_variant_context(self):
        print("aap")
