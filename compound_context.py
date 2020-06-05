"""Object to store linked context."""


class CompoundContext:
    """Stores information about linked variant contexts."""

    def __init__(self, first_varcon, second_varcon):
        self.variant_contexts = [first_varcon, second_varcon]
        self.priority_level = 64
        self.priority_label = ""

    def get_variant_contexts(self):
        """Return the variant contexts forming a compound context.

        Returns
        -------
        self.variant_contexts : list of VariantContext
            Variant contexts forming a
        """
        return self.variant_contexts

    def get_num_of_variant_contexts(self):
        """Return the number of variant contexts forming a compound context.

        Returns
        -------
        int
            The number of variant contexts forming a compound context
        """
        return len(self.variant_contexts)

    def add_variant_context(self, varcon_toadd):
        """Add a provided VariantContext to the compound context.

        Parameters
        ----------
        varcon_toadd : VariantContext
            self.variant_contexts.append(varcon_toadd)
        """
        if varcon_toadd is not None:
            self.variant_contexts.append(varcon_toadd)

    def add_variant_contexts(self, varcons_to_add):
        """Add a list of variant contexts to a compound context.

        Parameters
        ----------
        varcons_to_add : list of VariantContext
            List of variant contexts to add to the compound
        """
        if varcons_to_add is not None:
            varcons_to_add = [varcon for varcon in varcons_to_add if varcon]
            if len(varcons_to_add) >= 1:
                self.variant_contexts.extend(varcons_to_add)

    # TODO:
# =============================================================================
#     def determine_priority_level(self):
#         """WIP."""
#         print("test")
#
#     def remove_variant_context(self):
#         """WIP."""
#         print("test")
#
#
#     def get_variant_context_index(self, variant_context):
#         """Determine and return the index in the list of variant contexts.
#
#         This can be used to remove the proper variant context from the compound context.
#
#         Parameters
#         ----------
#         variant_context : VariantContext
#
#         Returns
#         -------
#         int
#             Index position in list of variant contexts
#         """
#         print("test")
# =============================================================================
