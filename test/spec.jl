using Test
import CryptoGroups: Specs
import Specs: gn_basis_exist, gn_basis_representation_rule



@test gn_basis_exist(4, 3) == true

# C.1 Table of GNB in X9.62
@test gn_basis_representation_rule(161) == 6
@test gn_basis_representation_rule(185) == 8
@test gn_basis_representation_rule(186) == 2
@test gn_basis_representation_rule(190) == 10
@test gn_basis_representation_rule(300) == 19
@test gn_basis_representation_rule(487) == 4
@test gn_basis_representation_rule(628) == 7
@test gn_basis_representation_rule(1380) == 1
@test gn_basis_representation_rule(1703) == 2

