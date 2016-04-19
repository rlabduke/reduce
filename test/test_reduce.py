#!/usr/bin/python

###############################################################################
#
#  This script tests various functions of reduce. I encourage you to add tests
#  as needs arise. When doing so, please update the info on what is ran here
#  below. The script uses python's unittest framework to compare reduce output
#  pdb strings to expected pdb strings which are stored as strings in this
#  script.
#
#  The following is a list and description of the current tests ran herein:
#    - test_default_reduce       tests default reduce
#                                Compares output to default_answer
#    - test_flip_reduce          test reduce's flip functionality using the
#                                -flip flag,
#                                Compares output to flip_answer
#    - test_nuc_reduce           tests the addition of nuclear-centered Hs.
#                                Compares output to nuc_answer
#
# TODO : Add N and H to the flip test
#
###############################################################################

import unittest
import os
import subprocess

hs = 'you must go in the tst script and change the %s variable to point to %s'
reduce_exe = '../reduce_src/reduce'
assert os.path.exists(reduce_exe), hs % ('reduce_exe','a reduce executable')
het_dict_path = '../reduce_het_dict.txt'
assert os.path.exists(het_dict_path),hs%('het_dict_path','reduce_het_dict.txt')

__author__ = 'bhintze'

default_pdb_str = '''\
CRYST1   42.360   54.790  111.990  90.00  90.00  90.00 P 21 21 21    4
SCALE1      0.023607  0.000000  0.000000        0.00000
SCALE2      0.000000  0.018252  0.000000        0.00000
ATOM    312  N   THR A  68      33.568  21.278  46.864  1.00 44.97           N
ATOM    313  CA  THR A  68      32.890  20.190  47.567  1.00 46.61           C
ATOM    314  C   THR A  68      32.844  20.354  49.087  1.00 48.10           C
ATOM    315  O   THR A  68      32.371  19.466  49.797  1.00 50.47           O
ATOM    316  CB  THR A  68      33.558  18.847  47.259  1.00 46.45           C
ATOM    317  OG1 THR A  68      34.862  18.831  47.845  1.00 54.25           O
ATOM    318  CG2 THR A  68      33.679  18.648  45.760  1.00 45.17           C
ATOM    319  N   GLY A  69      33.330  21.483  49.591  1.00 46.02           N
ATOM    320  CA  GLY A  69      33.272  21.741  51.020  1.00 45.56           C
ATOM    321  C   GLY A  69      31.908  22.257  51.441  1.00 49.19           C
ATOM    322  O   GLY A  69      30.977  22.314  50.636  1.00 47.01           O
ATOM    323  N   GLY A  70      31.779  22.626  52.712  1.00 50.78           N
ATOM    324  CA  GLY A  70      30.561  23.259  53.190  1.00 47.53           C
ATOM    325  C   GLY A  70      29.348  22.361  53.314  1.00 51.12           C
ATOM    326  O   GLY A  70      28.258  22.831  53.644  1.00 49.03           O
ATOM    327  N   ASN A  71      29.534  21.069  53.061  1.00 58.93           N
ATOM    328  CA  ASN A  71      28.417  20.129  53.040  1.00 65.14           C
ATOM    329  C   ASN A  71      27.932  19.713  54.427  1.00 74.93           C
ATOM    330  O   ASN A  71      28.629  19.904  55.427  1.00 77.88           O
ATOM    331  CB  ASN A  71      28.791  18.894  52.225  1.00 61.19           C
ATOM    332  CG  ASN A  71      28.608  19.110  50.735  1.00 58.59           C
ATOM    333  OD1 ASN A  71      27.771  18.468  50.099  1.00 62.04           O
ATOM    334  ND2 ASN A  71      29.382  20.033  50.174  1.00 58.77           N
ATOM    335  N   ILE A  72      26.737  19.123  54.447  1.00 76.95           N
ATOM    336  CA  ILE A  72      25.929  18.894  55.651  1.00 80.49           C
ATOM    337  C   ILE A  72      26.679  18.458  56.914  1.00 91.67           C
ATOM    338  O   ILE A  72      27.103  19.298  57.714  1.00 90.79           O
ATOM    339  N   ARG A  73      26.813  17.144  57.090  1.00 96.55           N
ATOM    340  CA  ARG A  73      27.329  16.556  58.318  1.00 95.73           C
ATOM    341  C   ARG A  73      28.678  17.073  58.781  1.00 97.01           C
ATOM    342  O   ARG A  73      28.915  17.221  59.983  1.00 95.84           O
ATOM    343  N   ASN A  74      29.564  17.343  57.826  1.00 93.35           N
ATOM    344  CA  ASN A  74      30.868  17.899  58.135  1.00 93.98           C
ATOM    345  C   ASN A  74      30.791  19.397  58.361  1.00 90.89           C
ATOM    346  O   ASN A  74      31.261  20.183  57.536  1.00 80.51           O
ATOM    347  N   ASP A  75      30.194  19.790  59.484  1.00 91.06           N
ATOM    348  CA  ASP A  75      29.985  21.193  59.783  1.00 82.35           C
ATOM    349  C   ASP A  75      29.897  21.516  61.260  1.00 79.45           C
ATOM    350  O   ASP A  75      29.021  22.271  61.686  1.00 75.34           O
ATOM    351  N   ASP A  76      30.806  20.947  62.046  1.00 82.28           N
ATOM    352  CA  ASP A  76      30.916  21.299  63.449  1.00 79.21           C
ATOM    353  C   ASP A  76      31.526  22.681  63.584  1.00 76.33           C
ATOM    354  O   ASP A  76      31.412  23.327  64.622  1.00 79.14           O
ATOM    355  N   LYS A  77      32.175  23.136  62.515  1.00 77.70           N
ATOM    356  CA  LYS A  77      32.808  24.442  62.495  1.00 74.51           C
ATOM    357  C   LYS A  77      32.161  25.410  61.518  1.00 68.95           C
ATOM    358  O   LYS A  77      32.764  26.412  61.132  1.00 62.03           O
ATOM    359  N   TYR A  78      30.936  25.104  61.103  1.00 64.80           N
ATOM    360  CA  TYR A  78      30.137  26.048  60.333  1.00 56.64           C
ATOM    361  C   TYR A  78      29.133  26.715  61.259  1.00 48.87           C
ATOM    362  O   TYR A  78      28.369  27.587  60.844  1.00 44.48           O
ATOM    363  CB  TYR A  78      29.417  25.352  59.173  1.00 57.64           C
ATOM    364  CG  TYR A  78      30.352  24.788  58.130  1.00 61.82           C
ATOM    365  CD1 TYR A  78      31.582  25.388  57.872  1.00 63.84           C
ATOM    366  CD2 TYR A  78      30.012  23.654  57.406  1.00 64.20           C
ATOM    367  CE1 TYR A  78      32.446  24.873  56.916  1.00 58.36           C
ATOM    368  CE2 TYR A  78      30.869  23.128  56.458  1.00 65.24           C
ATOM    369  CZ  TYR A  78      32.083  23.743  56.212  1.00 61.98           C
ATOM    370  OH  TYR A  78      32.941  23.224  55.265  1.00 63.38           O
ATOM    371  N   THR A  79      29.156  26.301  62.523  1.00 49.22           N
ATOM    372  CA  THR A  79      28.229  26.804  63.524  1.00 42.21           C
ATOM    373  C   THR A  79      28.321  28.320  63.691  1.00 38.81           C
ATOM    374  O   THR A  79      27.364  28.958  64.127  1.00 40.18           O
ATOM    375  CB  THR A  79      28.468  26.130  64.895  1.00 48.93           C
ATOM    376  OG1 THR A  79      29.771  26.471  65.381  1.00 55.67           O
ATOM    377  CG2 THR A  79      28.367  24.620  64.768  1.00 56.67           C
ATOM    378  N   HIS A  80      29.463  28.905  63.347  1.00 33.78           N
ATOM    379  CA  HIS A  80      29.612  30.349  63.476  1.00 31.59           C
ATOM    380  C   HIS A  80      28.915  31.110  62.336  1.00 28.85           C
ATOM    381  O   HIS A  80      28.790  32.333  62.393  1.00 28.35           O
ATOM    382  CB  HIS A  80      31.095  30.737  63.556  1.00 36.80           C
ATOM    383  CG  HIS A  80      31.885  30.403  62.328  1.00 44.69           C
ATOM    384  ND1 HIS A  80      31.515  29.411  61.446  1.00 46.50           N
ATOM    385  CD2 HIS A  80      33.036  30.929  61.843  1.00 44.38           C
ATOM    386  CE1 HIS A  80      32.398  29.346  60.463  1.00 43.56           C
ATOM    387  NE2 HIS A  80      33.329  30.257  60.679  1.00 44.44           N
ATOM    388  N   PHE A  81      28.455  30.392  61.314  1.00 26.91           N
ATOM    389  CA  PHE A  81      27.611  31.013  60.293  1.00 26.80           C
ATOM    390  C   PHE A  81      26.207  31.294  60.819  1.00 27.58           C
ATOM    391  O   PHE A  81      25.485  32.115  60.260  1.00 23.87           O
ATOM    392  CB  PHE A  81      27.512  30.132  59.040  1.00 26.61           C
ATOM    393  CG  PHE A  81      28.820  29.946  58.313  1.00 25.55           C
ATOM    394  CD1 PHE A  81      29.045  28.809  57.565  1.00 27.17           C
ATOM    395  CD2 PHE A  81      29.815  30.911  58.372  1.00 26.18           C
ATOM    396  CE1 PHE A  81      30.242  28.624  56.897  1.00 27.35           C
ATOM    397  CE2 PHE A  81      31.018  30.732  57.703  1.00 24.11           C
ATOM    398  CZ  PHE A  81      31.229  29.587  56.967  1.00 26.35           C
ATOM    399  N   PHE A  82      25.818  30.632  61.909  1.00 26.96           N
ATOM    400  CA  PHE A  82      24.448  30.780  62.397  1.00 28.96           C
ATOM    401  C   PHE A  82      24.398  31.122  63.878  1.00 29.75           C
ATOM    402  O   PHE A  82      25.280  30.748  64.630  1.00 34.14           O
ATOM    403  CB  PHE A  82      23.654  29.503  62.128  1.00 32.61           C
ATOM    404  CG  PHE A  82      23.702  29.060  60.693  1.00 32.62           C
ATOM    405  CD1 PHE A  82      24.535  28.027  60.297  1.00 31.95           C
ATOM    406  CD2 PHE A  82      22.930  29.691  59.738  1.00 28.66           C
ATOM    407  CE1 PHE A  82      24.584  27.632  58.973  1.00 30.64           C
ATOM    408  CE2 PHE A  82      22.978  29.296  58.407  1.00 30.13           C
ATOM    409  CZ  PHE A  82      23.805  28.268  58.032  1.00 30.16           C
ATOM    410  N   SER A  83      23.367  31.856  64.279  1.00 28.83           N
ATOM    411  CA  SER A  83      23.146  32.159  65.686  1.00 33.49           C
ATOM    412  C   SER A  83      22.802  30.880  66.450  1.00 39.39           C
ATOM    413  O   SER A  83      22.535  29.841  65.854  1.00 37.63           O
ATOM    414  CB  SER A  83      22.026  33.183  65.849  1.00 32.43           C
ATOM    415  OG  SER A  83      20.768  32.557  65.653  1.00 44.84           O
ATOM    416  N   GLY A  84      22.812  30.954  67.775  1.00 46.57           N
ATOM    417  CA  GLY A  84      22.555  29.775  68.586  1.00 49.01           C
ATOM    418  C   GLY A  84      21.110  29.319  68.510  1.00 46.09           C
ATOM    419  O   GLY A  84      20.830  28.124  68.567  1.00 53.68           O
HETATM 1739  C1  EOH A 717      34.026  19.250  53.995  0.90 63.57           C
HETATM 1740  C2  EOH A 717      33.701  20.048  55.258  0.90 61.11           C
HETATM 1741  O   EOH A 717      32.822  18.995  53.269  0.90 65.41           O
HETATM 1755  C   MOH A 724      24.486  26.604  64.638  1.00 43.73           C
HETATM 1756  O   MOH A 724      25.220  27.577  65.381  1.00 45.19           O
END
'''

default_answer = '''\
USER  MOD reduce.3.24.130724 H: found=0, std=0, add=78, rem=0, adj=5
CRYST1   42.360   54.790  111.990  90.00  90.00  90.00 P 21 21 21    4
SCALE1      0.023607  0.000000  0.000000        0.00000
SCALE2      0.000000  0.018252  0.000000        0.00000
USER  MOD -----------------------------------------------------------------
USER  MOD scores for adjustable sidechains, with "set" totals for H,N and Q
USER  MOD "o" means original, "f" means flipped, "180deg" is methyl default
USER  MOD "!" flags a clash with an overlap of 0.40A or greater
USER  MOD flip categories: "K"=keep, "C"=clashes, "X"=uncertain, "F"=flip
USER  MOD Single : A  68 THR OG1 :   rot  -40:sc= 0.00645
USER  MOD Single : A  78 TYR OH  :   rot   63:sc=  -0.281
USER  MOD Single : A  79 THR OG1 :   rot  180:sc=       0
USER  MOD Single : A  83 SER OG  :   rot  180:sc=       0
USER  MOD Single : A 717 EOH O   :   rot  180:sc=       0
USER  MOD -----------------------------------------------------------------
ATOM    312  N   THR A  68      33.568  21.278  46.864  1.00 44.97           N
ATOM    313  CA  THR A  68      32.890  20.190  47.567  1.00 46.61           C
ATOM    314  C   THR A  68      32.844  20.354  49.087  1.00 48.10           C
ATOM    315  O   THR A  68      32.371  19.466  49.797  1.00 50.47           O
ATOM    316  CB  THR A  68      33.558  18.847  47.259  1.00 46.45           C
ATOM    317  OG1 THR A  68      34.862  18.831  47.845  1.00 54.25           O
ATOM    318  CG2 THR A  68      33.679  18.648  45.760  1.00 45.17           C
ATOM      0  HA  THR A  68      31.977  20.217  47.240  1.00 46.61           H   new
ATOM      0  HB  THR A  68      33.017  18.131  47.627  1.00 46.45           H   new
ATOM      0  HG1 THR A  68      35.224  19.582  47.745  1.00 54.25           H   new
ATOM      0 HG21 THR A  68      34.103  17.795  45.580  1.00 45.17           H   new
ATOM      0 HG22 THR A  68      32.796  18.660  45.359  1.00 45.17           H   new
ATOM      0 HG23 THR A  68      34.216  19.362  45.381  1.00 45.17           H   new
ATOM    319  N   GLY A  69      33.330  21.483  49.591  1.00 46.02           N
ATOM    320  CA  GLY A  69      33.272  21.741  51.020  1.00 45.56           C
ATOM    321  C   GLY A  69      31.908  22.257  51.441  1.00 49.19           C
ATOM    322  O   GLY A  69      30.977  22.314  50.636  1.00 47.01           O
ATOM      0  H   GLY A  69      33.694  22.107  49.125  1.00 46.02           H   new
ATOM      0  HA2 GLY A  69      33.474  20.926  51.505  1.00 45.56           H   new
ATOM      0  HA3 GLY A  69      33.952  22.390  51.260  1.00 45.56           H   new
ATOM    323  N   GLY A  70      31.779  22.626  52.712  1.00 50.78           N
ATOM    324  CA  GLY A  70      30.561  23.259  53.190  1.00 47.53           C
ATOM    325  C   GLY A  70      29.348  22.361  53.314  1.00 51.12           C
ATOM    326  O   GLY A  70      28.258  22.831  53.644  1.00 49.03           O
ATOM      0  H   GLY A  70      32.386  22.518  53.311  1.00 50.78           H   new
ATOM      0  HA2 GLY A  70      30.741  23.650  54.059  1.00 47.53           H   new
ATOM      0  HA3 GLY A  70      30.341  23.989  52.590  1.00 47.53           H   new
ATOM    327  N   ASN A  71      29.534  21.069  53.061  1.00 58.93           N
ATOM    328  CA  ASN A  71      28.417  20.129  53.040  1.00 65.14           C
ATOM    329  C   ASN A  71      27.932  19.713  54.427  1.00 74.93           C
ATOM    330  O   ASN A  71      28.629  19.904  55.427  1.00 77.88           O
ATOM    331  CB  ASN A  71      28.791  18.894  52.225  1.00 61.19           C
ATOM    332  CG  ASN A  71      28.608  19.110  50.735  1.00 58.59           C
ATOM    333  OD1 ASN A  71      27.771  18.468  50.099  1.00 62.04           O
ATOM    334  ND2 ASN A  71      29.382  20.033  50.174  1.00 58.77           N
ATOM      0  H   ASN A  71      30.301  20.716  52.899  1.00 58.93           H   new
ATOM      0  HA  ASN A  71      27.676  20.597  52.624  1.00 65.14           H   new
ATOM      0  HB2 ASN A  71      29.715  18.658  52.404  1.00 61.19           H   new
ATOM      0  HB3 ASN A  71      28.246  18.144  52.510  1.00 61.19           H   new
ATOM      0 HD21 ASN A  71      29.309  20.202  49.334  1.00 58.77           H   new
ATOM      0 HD22 ASN A  71      29.955  20.461  50.651  1.00 58.77           H   new
ATOM    335  N   ILE A  72      26.737  19.123  54.447  1.00 76.95           N
ATOM    336  CA  ILE A  72      25.929  18.894  55.651  1.00 80.49           C
ATOM    337  C   ILE A  72      26.679  18.458  56.914  1.00 91.67           C
ATOM    338  O   ILE A  72      27.103  19.298  57.714  1.00 90.79           O
ATOM      0  H   ILE A  72      26.358  18.833  53.732  1.00 76.95           H   new
ATOM    339  N   ARG A  73      26.813  17.144  57.090  1.00 96.55           N
ATOM    340  CA  ARG A  73      27.329  16.556  58.318  1.00 95.73           C
ATOM    341  C   ARG A  73      28.678  17.073  58.781  1.00 97.01           C
ATOM    342  O   ARG A  73      28.915  17.221  59.983  1.00 95.84           O
ATOM      0  H   ARG A  73      26.603  16.564  56.490  1.00 96.55           H   new
ATOM    343  N   ASN A  74      29.564  17.343  57.826  1.00 93.35           N
ATOM    344  CA  ASN A  74      30.868  17.899  58.135  1.00 93.98           C
ATOM    345  C   ASN A  74      30.791  19.397  58.361  1.00 90.89           C
ATOM    346  O   ASN A  74      31.261  20.183  57.536  1.00 80.51           O
ATOM      0  H   ASN A  74      29.424  17.209  56.988  1.00 93.35           H   new
ATOM    347  N   ASP A  75      30.194  19.790  59.484  1.00 91.06           N
ATOM    348  CA  ASP A  75      29.985  21.193  59.783  1.00 82.35           C
ATOM    349  C   ASP A  75      29.897  21.516  61.260  1.00 79.45           C
ATOM    350  O   ASP A  75      29.021  22.271  61.686  1.00 75.34           O
ATOM      0  H   ASP A  75      29.902  19.251  60.087  1.00 91.06           H   new
ATOM    351  N   ASP A  76      30.806  20.947  62.046  1.00 82.28           N
ATOM    352  CA  ASP A  76      30.916  21.299  63.449  1.00 79.21           C
ATOM    353  C   ASP A  76      31.526  22.681  63.584  1.00 76.33           C
ATOM    354  O   ASP A  76      31.412  23.327  64.622  1.00 79.14           O
ATOM      0  H   ASP A  76      31.369  20.353  61.781  1.00 82.28           H   new
ATOM    355  N   LYS A  77      32.175  23.136  62.515  1.00 77.70           N
ATOM    356  CA  LYS A  77      32.808  24.442  62.495  1.00 74.51           C
ATOM    357  C   LYS A  77      32.161  25.410  61.518  1.00 68.95           C
ATOM    358  O   LYS A  77      32.764  26.412  61.132  1.00 62.03           O
ATOM      0  H   LYS A  77      32.258  22.693  61.783  1.00 77.70           H   new
ATOM    359  N   TYR A  78      30.936  25.104  61.103  1.00 64.80           N
ATOM    360  CA  TYR A  78      30.137  26.048  60.333  1.00 56.64           C
ATOM    361  C   TYR A  78      29.133  26.715  61.259  1.00 48.87           C
ATOM    362  O   TYR A  78      28.369  27.587  60.844  1.00 44.48           O
ATOM    363  CB  TYR A  78      29.417  25.352  59.173  1.00 57.64           C
ATOM    364  CG  TYR A  78      30.352  24.788  58.130  1.00 61.82           C
ATOM    365  CD1 TYR A  78      31.582  25.388  57.872  1.00 63.84           C
ATOM    366  CD2 TYR A  78      30.012  23.654  57.406  1.00 64.20           C
ATOM    367  CE1 TYR A  78      32.446  24.873  56.916  1.00 58.36           C
ATOM    368  CE2 TYR A  78      30.869  23.128  56.458  1.00 65.24           C
ATOM    369  CZ  TYR A  78      32.083  23.743  56.212  1.00 61.98           C
ATOM    370  OH  TYR A  78      32.941  23.224  55.265  1.00 63.38           O
ATOM      0  H   TYR A  78      30.549  24.352  61.258  1.00 64.80           H   new
ATOM      0  HA  TYR A  78      30.724  26.717  59.949  1.00 56.64           H   new
ATOM      0  HB2 TYR A  78      28.869  24.634  59.527  1.00 57.64           H   new
ATOM      0  HB3 TYR A  78      28.816  25.985  58.749  1.00 57.64           H   new
ATOM      0  HD1 TYR A  78      31.829  26.147  58.349  1.00 63.84           H   new
ATOM      0  HD2 TYR A  78      29.193  23.241  57.561  1.00 64.20           H   new
ATOM      0  HE1 TYR A  78      33.263  25.286  56.751  1.00 58.36           H   new
ATOM      0  HE2 TYR A  78      30.630  22.363  55.987  1.00 65.24           H   new
ATOM      0  HH  TYR A  78      33.653  22.977  55.637  1.00 63.38           H   new
ATOM    371  N   THR A  79      29.156  26.301  62.523  1.00 49.22           N
ATOM    372  CA  THR A  79      28.229  26.804  63.524  1.00 42.21           C
ATOM    373  C   THR A  79      28.321  28.320  63.691  1.00 38.81           C
ATOM    374  O   THR A  79      27.364  28.958  64.127  1.00 40.18           O
ATOM    375  CB  THR A  79      28.468  26.130  64.895  1.00 48.93           C
ATOM    376  OG1 THR A  79      29.771  26.471  65.381  1.00 55.67           O
ATOM    377  CG2 THR A  79      28.367  24.620  64.768  1.00 56.67           C
ATOM      0  H   THR A  79      29.713  25.718  62.822  1.00 49.22           H   new
ATOM      0  HA  THR A  79      27.340  26.585  63.203  1.00 42.21           H   new
ATOM      0  HB  THR A  79      27.791  26.445  65.515  1.00 48.93           H   new
ATOM      0  HG1 THR A  79      29.897  26.105  66.126  1.00 55.67           H   new
ATOM      0 HG21 THR A  79      28.519  24.211  65.634  1.00 56.67           H   new
ATOM      0 HG22 THR A  79      27.483  24.380  64.448  1.00 56.67           H   new
ATOM      0 HG23 THR A  79      29.034  24.303  64.140  1.00 56.67           H   new
ATOM    378  N   HIS A  80      29.463  28.905  63.347  1.00 33.78           N
ATOM    379  CA  HIS A  80      29.612  30.349  63.476  1.00 31.59           C
ATOM    380  C   HIS A  80      28.915  31.110  62.336  1.00 28.85           C
ATOM    381  O   HIS A  80      28.790  32.333  62.393  1.00 28.35           O
ATOM    382  CB  HIS A  80      31.095  30.737  63.556  1.00 36.80           C
ATOM    383  CG  HIS A  80      31.885  30.403  62.328  1.00 44.69           C
ATOM    384  ND1 HIS A  80      31.515  29.411  61.446  1.00 46.50           N
ATOM    385  CD2 HIS A  80      33.036  30.929  61.843  1.00 44.38           C
ATOM    386  CE1 HIS A  80      32.398  29.346  60.463  1.00 43.56           C
ATOM    387  NE2 HIS A  80      33.329  30.257  60.679  1.00 44.44           N
ATOM      0  H   HIS A  80      30.153  28.492  63.042  1.00 33.78           H   new
ATOM      0  HA  HIS A  80      29.176  30.607  64.303  1.00 31.59           H   new
ATOM      0  HB2 HIS A  80      31.161  31.691  63.721  1.00 36.80           H   new
ATOM      0  HB3 HIS A  80      31.495  30.289  64.318  1.00 36.80           H   new
ATOM      0  HD2 HIS A  80      33.534  31.615  62.225  1.00 44.38           H   new
ATOM      0  HE1 HIS A  80      32.368  28.759  59.743  1.00 43.56           H   new
ATOM    388  N   PHE A  81      28.455  30.392  61.314  1.00 26.91           N
ATOM    389  CA  PHE A  81      27.611  31.013  60.293  1.00 26.80           C
ATOM    390  C   PHE A  81      26.207  31.294  60.819  1.00 27.58           C
ATOM    391  O   PHE A  81      25.485  32.115  60.260  1.00 23.87           O
ATOM    392  CB  PHE A  81      27.512  30.132  59.040  1.00 26.61           C
ATOM    393  CG  PHE A  81      28.820  29.946  58.313  1.00 25.55           C
ATOM    394  CD1 PHE A  81      29.045  28.809  57.565  1.00 27.17           C
ATOM    395  CD2 PHE A  81      29.815  30.911  58.372  1.00 26.18           C
ATOM    396  CE1 PHE A  81      30.242  28.624  56.897  1.00 27.35           C
ATOM    397  CE2 PHE A  81      31.018  30.732  57.703  1.00 24.11           C
ATOM    398  CZ  PHE A  81      31.229  29.587  56.967  1.00 26.35           C
ATOM      0  H   PHE A  81      28.616  29.556  61.193  1.00 26.91           H   new
ATOM      0  HA  PHE A  81      28.033  31.855  60.059  1.00 26.80           H   new
ATOM      0  HB2 PHE A  81      27.168  29.262  59.295  1.00 26.61           H   new
ATOM      0  HB3 PHE A  81      26.867  30.524  58.431  1.00 26.61           H   new
ATOM      0  HD1 PHE A  81      28.383  28.158  57.509  1.00 27.17           H   new
ATOM      0  HD2 PHE A  81      29.674  31.687  58.865  1.00 26.18           H   new
ATOM      0  HE1 PHE A  81      30.382  27.850  56.400  1.00 27.35           H   new
ATOM      0  HE2 PHE A  81      31.680  31.384  57.752  1.00 24.11           H   new
ATOM      0  HZ  PHE A  81      32.034  29.463  56.518  1.00 26.35           H   new
ATOM    399  N   PHE A  82      25.818  30.632  61.909  1.00 26.96           N
ATOM    400  CA  PHE A  82      24.448  30.780  62.397  1.00 28.96           C
ATOM    401  C   PHE A  82      24.398  31.122  63.878  1.00 29.75           C
ATOM    402  O   PHE A  82      25.280  30.748  64.630  1.00 34.14           O
ATOM    403  CB  PHE A  82      23.654  29.503  62.128  1.00 32.61           C
ATOM    404  CG  PHE A  82      23.702  29.060  60.693  1.00 32.62           C
ATOM    405  CD1 PHE A  82      24.535  28.027  60.297  1.00 31.95           C
ATOM    406  CD2 PHE A  82      22.930  29.691  59.738  1.00 28.66           C
ATOM    407  CE1 PHE A  82      24.584  27.632  58.973  1.00 30.64           C
ATOM    408  CE2 PHE A  82      22.978  29.296  58.407  1.00 30.13           C
ATOM    409  CZ  PHE A  82      23.805  28.268  58.032  1.00 30.16           C
ATOM      0  H   PHE A  82      26.318  30.105  62.370  1.00 26.96           H   new
ATOM      0  HA  PHE A  82      24.048  31.521  61.915  1.00 28.96           H   new
ATOM      0  HB2 PHE A  82      23.998  28.792  62.691  1.00 32.61           H   new
ATOM      0  HB3 PHE A  82      22.729  29.645  62.385  1.00 32.61           H   new
ATOM      0  HD1 PHE A  82      25.066  27.595  60.927  1.00 31.95           H   new
ATOM      0  HD2 PHE A  82      22.370  30.390  59.989  1.00 28.66           H   new
ATOM      0  HE1 PHE A  82      25.144  26.935  58.717  1.00 30.64           H   new
ATOM      0  HE2 PHE A  82      22.452  29.727  57.773  1.00 30.13           H   new
ATOM      0  HZ  PHE A  82      23.840  28.000  57.142  1.00 30.16           H   new
ATOM    410  N   SER A  83      23.367  31.856  64.279  1.00 28.83           N
ATOM    411  CA  SER A  83      23.146  32.159  65.686  1.00 33.49           C
ATOM    412  C   SER A  83      22.802  30.880  66.450  1.00 39.39           C
ATOM    413  O   SER A  83      22.535  29.841  65.854  1.00 37.63           O
ATOM    414  CB  SER A  83      22.026  33.183  65.849  1.00 32.43           C
ATOM    415  OG  SER A  83      20.768  32.557  65.653  1.00 44.84           O
ATOM      0  H   SER A  83      22.780  32.190  63.747  1.00 28.83           H   new
ATOM      0  HA  SER A  83      23.962  32.537  66.050  1.00 33.49           H   new
ATOM      0  HB2 SER A  83      22.065  33.579  66.734  1.00 32.43           H   new
ATOM      0  HB3 SER A  83      22.141  33.904  65.210  1.00 32.43           H   new
ATOM      0  HG  SER A  83      20.155  33.124  65.745  1.00 44.84           H   new
ATOM    416  N   GLY A  84      22.812  30.954  67.775  1.00 46.57           N
ATOM    417  CA  GLY A  84      22.555  29.775  68.586  1.00 49.01           C
ATOM    418  C   GLY A  84      21.110  29.319  68.510  1.00 46.09           C
ATOM    419  O   GLY A  84      20.830  28.124  68.567  1.00 53.68           O
ATOM      0  H   GLY A  84      22.965  31.673  68.221  1.00 46.57           H   new
ATOM      0  HA2 GLY A  84      23.134  29.054  68.295  1.00 49.01           H   new
ATOM      0  HA3 GLY A  84      22.782  29.966  69.509  1.00 49.01           H   new
HETATM 1739  C1  EOH A 717      34.026  19.250  53.995  0.90 63.57           C
HETATM 1740  C2  EOH A 717      33.701  20.048  55.258  0.90 61.11           C
HETATM 1741  O   EOH A 717      32.822  18.995  53.269  0.90 65.41           O
HETATM    0 3H2  EOH A 717      33.288  20.890  55.011  0.90 61.11           H   new
HETATM    0 2H2  EOH A 717      33.090  19.539  55.814  0.90 61.11           H   new
HETATM    0 2H1  EOH A 717      34.651  19.743  53.441  0.90 63.57           H   new
HETATM    0 1H2  EOH A 717      34.518  20.221  55.751  0.90 61.11           H   new
HETATM    0 1H1  EOH A 717      34.456  18.413  54.231  0.90 63.57           H   new
HETATM    0  HO  EOH A 717      33.001  18.558  52.574  0.90 65.41           H   new
HETATM 1755  C   MOH A 724      24.486  26.604  64.638  1.00 43.73           C
HETATM 1756  O   MOH A 724      25.220  27.577  65.381  1.00 45.19           O
END
'''

flip_test_str = '''\
ATOM    800  N   ASN B   3       6.076  89.982  59.411  1.00  4.63           N
ATOM    801  CA  ASN B   3       5.450  88.943  60.224  1.00  4.43           C
ATOM    802  C   ASN B   3       5.198  89.487  61.636  1.00  4.66           C
ATOM    803  O   ASN B   3       4.110  89.265  62.214  1.00  4.93           O
ATOM    804  CB  ASN B   3       6.301  87.657  60.231  1.00  4.34           C
ATOM    805  CG  ASN B   3       5.552  86.456  60.791  1.00  3.94           C
ATOM    806  OD1 ASN B   3       4.451  86.590  61.345  1.00  2.54           O
ATOM    807  ND2 ASN B   3       6.148  85.267  60.647  1.00  4.44           N
ATOM   1837  N   HIS B   4       6.171  90.222  62.163  1.00  2.00           N
ATOM   1838  CA  HIS B   4       6.033  90.827  63.489  1.00  2.00           C
ATOM   1839  C   HIS B   4       4.816  91.761  63.577  1.00  2.00           C
ATOM   1840  O   HIS B   4       4.170  91.880  64.649  1.00  2.00           O
ATOM   1841  CB  HIS B   4       7.302  91.606  63.860  1.00  2.00           C
ATOM   1842  CG  HIS B   4       8.501  90.745  64.094  1.00  2.00           C
ATOM   1843  ND1 HIS B   4       9.700  90.948  63.443  1.00  2.00           N
ATOM   1844  CD2 HIS B   4       8.696  89.694  64.924  1.00  2.00           C
ATOM   1845  CE1 HIS B   4      10.576  90.049  63.853  1.00  2.00           C
ATOM   1846  NE2 HIS B   4       9.995  89.282  64.757  1.00  2.00           N
ATOM    151  N   ARG B   5       4.459  92.414  62.487  1.00 18.73           N
ATOM    152  CA  ARG B   5       3.320  93.335  62.497  1.00 19.66           C
ATOM    153  C   ARG B   5       2.032  92.680  62.014  1.00 19.66           C
ATOM    154  O   ARG B   5       0.950  93.273  62.111  1.00 21.95           O
ATOM    155  CB  ARG B   5       3.638  94.541  61.601  1.00 23.29           C
ATOM    156  CG  ARG B   5       4.881  95.354  62.034  1.00 23.24           C
ATOM    157  CD  ARG B   5       5.365  96.251  60.910  1.00 29.00           C
ATOM    158  NE  ARG B   5       6.449  97.202  61.226  1.00 32.94           N
ATOM    159  CZ  ARG B   5       7.728  96.894  61.473  1.00 34.16           C
ATOM    160  NH1 ARG B   5       8.157  95.628  61.480  1.00 30.53           N
ATOM    161  NH2 ARG B   5       8.606  97.884  61.648  1.00 35.22           N
ATOM    112  N   GLN B   9      10.186  91.359  59.521  1.00 15.79           N
ATOM    113  CA  GLN B   9       9.578  92.637  59.819  1.00 18.37           C
ATOM    114  C   GLN B   9       8.076  92.653  59.539  1.00 17.17           C
ATOM    115  O   GLN B   9       7.278  93.226  60.303  1.00 16.99           O
ATOM    116  CB  GLN B   9      10.255  93.750  59.036  1.00 19.44           C
ATOM    117  CG  GLN B   9      11.594  94.138  59.641  1.00 24.67           C
ATOM    118  CD  GLN B   9      11.403  94.696  61.052  0.00 27.67           C
ATOM    119  NE2 GLN B   9      10.724  95.708  61.264  1.00 31.74           O
ATOM    120  OE1 GLN B   9      11.983  94.024  62.020  1.00 31.80           N
HETATM 7573  O   HOH B 407       6.633  88.635  66.643  1.00 11.41           O
HETATM 7693  O   HOH B 527       3.491  85.743  61.660  1.00  4.12           O
'''

flip_answer = '''\
USER  MOD reduce.3.24.130724 H: found=0, std=0, add=33, rem=0, adj=3
USER  MOD -----------------------------------------------------------------
USER  MOD scores for adjustable sidechains, with "set" totals for H,N and Q
USER  MOD "o" means original, "f" means flipped, "180deg" is methyl default
USER  MOD "!" flags a clash with an overlap of 0.40A or greater
USER  MOD flip categories: "K"=keep, "C"=clashes, "X"=uncertain, "F"=flip
USER  MOD Single : B   3 ASN     :FLIP  amide:sc=   -4.98! C(o=-12!,f=-5!)
USER  MOD Single : B   4 HIS     :FLIP no HD1:sc=   0.634  F(o=-2.1!,f=0.63)
USER  MOD Single : B   9 GLN     :FLIP  amide:sc=   0.716  F(o=-4.8!,f=0.72)
USER  MOD -----------------------------------------------------------------
ATOM    800  N   ASN B   3       6.076  89.982  59.411  1.00  4.63           N
ATOM    801  CA  ASN B   3       5.450  88.943  60.224  1.00  4.43           C
ATOM    802  C   ASN B   3       5.198  89.487  61.636  1.00  4.66           C
ATOM    803  O   ASN B   3       4.110  89.265  62.214  1.00  4.93           O
ATOM    804  CB  ASN B   3       6.341  87.684  60.206  1.00  4.34           C
ATOM    805  CG  ASN B   3       5.632  86.451  60.747  1.00  3.94           C
ATOM    806  OD1 ASN B   3       6.138  85.324  60.640  1.00  2.54           O   flip
ATOM    807  ND2 ASN B   3       4.447  86.658  61.332  1.00  4.44           N   flip
ATOM      0  HA  ASN B   3       4.594  88.700  59.838  1.00  4.43           H   new
ATOM      0  HB2 ASN B   3       6.632  87.512  59.297  1.00  4.34           H   new
ATOM      0  HB3 ASN B   3       7.138  87.850  60.733  1.00  4.34           H   new
ATOM      0 HD21 ASN B   3       4.004  85.995  61.653  1.00  4.44           H   new
ATOM      0 HD22 ASN B   3       4.128  87.455  61.387  1.00  4.44           H   new
ATOM   1837  N   HIS B   4       6.171  90.222  62.163  1.00  2.00           N
ATOM   1838  CA  HIS B   4       6.033  90.827  63.489  1.00  2.00           C
ATOM   1839  C   HIS B   4       4.816  91.761  63.577  1.00  2.00           C
ATOM   1840  O   HIS B   4       4.170  91.880  64.649  1.00  2.00           O
ATOM   1841  CB  HIS B   4       7.265  91.675  63.832  1.00  2.00           C
ATOM   1842  CG  HIS B   4       8.531  90.890  63.957  1.00  2.00           C
ATOM   1843  ND1 HIS B   4       8.652  89.803  64.797  1.00  2.00           N   flip
ATOM   1844  CD2 HIS B   4       9.736  91.040  63.360  1.00  2.00           C   flip
ATOM   1845  CE1 HIS B   4       9.881  89.327  64.720  1.00  2.00           C   flip
ATOM   1846  NE2 HIS B   4      10.556  90.054  63.850  1.00  2.00           N   flip
ATOM      0  H   HIS B   4       6.920  90.384  61.772  1.00  2.00           H   new
ATOM      0  HA  HIS B   4       5.899  90.100  64.117  1.00  2.00           H   new
ATOM      0  HB2 HIS B   4       7.382  92.351  63.147  1.00  2.00           H   new
ATOM      0  HB3 HIS B   4       7.102  92.142  64.667  1.00  2.00           H   new
ATOM      0  HD2 HIS B   4       9.966  91.690  62.735  1.00  2.00           H   new
ATOM      0  HE1 HIS B   4      10.215  88.602  65.198  1.00  2.00           H   new
ATOM      0  HE2 HIS B   4      11.377  89.930  63.625  1.00  2.00           H   new
ATOM    151  N   ARG B   5       4.459  92.414  62.487  1.00 18.73           N
ATOM    152  CA  ARG B   5       3.320  93.335  62.497  1.00 19.66           C
ATOM    153  C   ARG B   5       2.032  92.680  62.014  1.00 19.66           C
ATOM    154  O   ARG B   5       0.950  93.273  62.111  1.00 21.95           O
ATOM    155  CB  ARG B   5       3.638  94.541  61.601  1.00 23.29           C
ATOM    156  CG  ARG B   5       4.881  95.354  62.034  1.00 23.24           C
ATOM    157  CD  ARG B   5       5.365  96.251  60.910  1.00 29.00           C
ATOM    158  NE  ARG B   5       6.449  97.202  61.226  1.00 32.94           N
ATOM    159  CZ  ARG B   5       7.728  96.894  61.473  1.00 34.16           C
ATOM    160  NH1 ARG B   5       8.157  95.628  61.480  1.00 30.53           N
ATOM    161  NH2 ARG B   5       8.606  97.884  61.648  1.00 35.22           N
ATOM      0  H   ARG B   5       4.857  92.344  61.728  1.00 18.73           H   new
ATOM      0  HA  ARG B   5       3.179  93.612  63.416  1.00 19.66           H   new
ATOM      0  HB2 ARG B   5       3.771  94.228  60.693  1.00 23.29           H   new
ATOM      0  HB3 ARG B   5       2.868  95.131  61.588  1.00 23.29           H   new
ATOM      0  HG2 ARG B   5       4.664  95.893  62.811  1.00 23.24           H   new
ATOM      0  HG3 ARG B   5       5.592  94.749  62.298  1.00 23.24           H   new
ATOM      0  HD2 ARG B   5       5.664  95.685  60.181  1.00 29.00           H   new
ATOM      0  HD3 ARG B   5       4.606  96.758  60.582  1.00 29.00           H   new
ATOM      0  HE  ARG B   5       6.238  98.035  61.254  1.00 32.94           H   new
ATOM      0 HH11 ARG B   5       7.606  94.986  61.323  1.00 30.53           H   new
ATOM      0 HH12 ARG B   5       8.984  95.455  61.641  1.00 30.53           H   new
ATOM      0 HH21 ARG B   5       8.346  98.702  61.601  1.00 35.22           H   new
ATOM      0 HH22 ARG B   5       9.432  97.703  61.808  1.00 35.22           H   new
ATOM    112  N   GLN B   9      10.186  91.359  59.521  1.00 15.79           N
ATOM    113  CA  GLN B   9       9.578  92.637  59.819  1.00 18.37           C
ATOM    114  C   GLN B   9       8.076  92.653  59.539  1.00 17.17           C
ATOM    115  O   GLN B   9       7.278  93.226  60.303  1.00 16.99           O
ATOM    116  CB  GLN B   9      10.229  93.787  59.068  1.00 19.44           C
ATOM    117  CG  GLN B   9      11.591  94.142  59.641  1.00 24.67           C
ATOM    118  CD  GLN B   9      11.455  94.628  61.084  0.00 27.67           C
ATOM    119  NE2 GLN B   9      11.988  94.028  62.025  1.00 31.74           O   flip
ATOM    120  OE1 GLN B   9      10.723  95.704  61.263  1.00 31.80           N   flip
ATOM      0  HA  GLN B   9       9.701  92.785  60.770  1.00 18.37           H   new
ATOM      0  HB2 GLN B   9      10.324  93.549  58.133  1.00 19.44           H   new
ATOM      0  HB3 GLN B   9       9.651  94.565  59.105  1.00 19.44           H   new
ATOM      0  HG2 GLN B   9      12.173  93.367  59.608  1.00 24.67           H   new
ATOM      0  HG3 GLN B   9      12.006  94.831  59.099  1.00 24.67           H   new
ATOM      0 HE21 GLN B   9      12.462  93.326  61.874  1.00 31.74           H   new
ATOM      0 HE22 GLN B   9      11.885  94.312  62.830  1.00 31.74           H   new
HETATM 7573  O   HOH B 407       6.633  88.635  66.643  1.00 11.41           O
HETATM 7693  O   HOH B 527       3.491  85.743  61.660  1.00  4.12           O
'''

nuc_answer =  '''\
USER  MOD reduce.3.24.130724 H: found=0, std=0, add=33, rem=0, adj=3
USER  MOD -----------------------------------------------------------------
USER  MOD scores for adjustable sidechains, with "set" totals for H,N and Q
USER  MOD "o" means original, "f" means flipped, "180deg" is methyl default
USER  MOD "!" flags a clash with an overlap of 0.40A or greater
USER  MOD flip categories: "K"=keep, "C"=clashes, "X"=uncertain, "F"=flip
USER  MOD Single : B   3 ASN     :FLIP  amide:sc=   -4.08! C(o=-13!,f=-4.1!)
USER  MOD Single : B   4 HIS     :FLIP no HE2:sc=   0.911  F(o=-2.8!,f=0.91)
USER  MOD Single : B   9 GLN     :FLIP  amide:sc=   0.688  F(o=-5.3!,f=0.69)
USER  MOD -----------------------------------------------------------------
ATOM    800  N   ASN B   3       6.076  89.982  59.411  1.00  4.63           N
ATOM    801  CA  ASN B   3       5.450  88.943  60.224  1.00  4.43           C
ATOM    802  C   ASN B   3       5.198  89.487  61.636  1.00  4.66           C
ATOM    803  O   ASN B   3       4.110  89.265  62.214  1.00  4.93           O
ATOM    804  CB  ASN B   3       6.341  87.684  60.206  1.00  4.34           C
ATOM    805  CG  ASN B   3       5.632  86.451  60.747  1.00  3.94           C
ATOM    806  OD1 ASN B   3       6.138  85.324  60.640  1.00  2.54           O   flip
ATOM    807  ND2 ASN B   3       4.447  86.658  61.332  1.00  4.44           N   flip
ATOM      0  HA  ASN B   3       4.488  88.669  59.790  1.00  4.43           H   new
ATOM      0  HB2 ASN B   3       6.668  87.491  59.184  1.00  4.34           H   new
ATOM      0  HB3 ASN B   3       7.237  87.871  60.798  1.00  4.34           H   new
ATOM      0 HD21 ASN B   3       3.921  85.871  61.713  1.00  4.44           H   new
ATOM      0 HD22 ASN B   3       4.069  87.603  61.397  1.00  4.44           H   new
ATOM   1837  N   HIS B   4       6.171  90.222  62.163  1.00  2.00           N
ATOM   1838  CA  HIS B   4       6.033  90.827  63.489  1.00  2.00           C
ATOM   1839  C   HIS B   4       4.816  91.761  63.577  1.00  2.00           C
ATOM   1840  O   HIS B   4       4.170  91.880  64.649  1.00  2.00           O
ATOM   1841  CB  HIS B   4       7.265  91.675  63.832  1.00  2.00           C
ATOM   1842  CG  HIS B   4       8.531  90.890  63.957  1.00  2.00           C
ATOM   1843  ND1 HIS B   4       8.652  89.803  64.797  1.00  2.00           N   flip
ATOM   1844  CD2 HIS B   4       9.736  91.040  63.360  1.00  2.00           C   flip
ATOM   1845  CE1 HIS B   4       9.881  89.327  64.720  1.00  2.00           C   flip
ATOM   1846  NE2 HIS B   4      10.556  90.054  63.850  1.00  2.00           N   flip
ATOM      0  H   HIS B   4       7.059  90.414  61.699  1.00  2.00           H   new
ATOM      0  HA  HIS B   4       5.883  90.010  64.195  1.00  2.00           H   new
ATOM      0  HB2 HIS B   4       7.396  92.434  63.062  1.00  2.00           H   new
ATOM      0  HB3 HIS B   4       7.082  92.200  64.770  1.00  2.00           H   new
ATOM      0  HD1 HIS B   4       7.909  89.426  65.386  1.00  2.00           H   new
ATOM      0  HD2 HIS B   4      10.003  91.795  62.635  1.00  2.00           H   new
ATOM      0  HE1 HIS B   4      10.268  88.485  65.275  1.00  2.00           H   new
ATOM    151  N   ARG B   5       4.459  92.414  62.487  1.00 18.73           N
ATOM    152  CA  ARG B   5       3.320  93.335  62.497  1.00 19.66           C
ATOM    153  C   ARG B   5       2.032  92.680  62.014  1.00 19.66           C
ATOM    154  O   ARG B   5       0.950  93.273  62.111  1.00 21.95           O
ATOM    155  CB  ARG B   5       3.638  94.541  61.601  1.00 23.29           C
ATOM    156  CG  ARG B   5       4.881  95.354  62.034  1.00 23.24           C
ATOM    157  CD  ARG B   5       5.365  96.251  60.910  1.00 29.00           C
ATOM    158  NE  ARG B   5       6.449  97.202  61.226  1.00 32.94           N
ATOM    159  CZ  ARG B   5       7.728  96.894  61.473  1.00 34.16           C
ATOM    160  NH1 ARG B   5       8.157  95.628  61.480  1.00 30.53           N
ATOM    161  NH2 ARG B   5       8.606  97.884  61.648  1.00 35.22           N
ATOM      0  H   ARG B   5       4.932  92.330  61.587  1.00 18.73           H   new
ATOM      0  HA  ARG B   5       3.162  93.646  63.530  1.00 19.66           H   new
ATOM      0  HB2 ARG B   5       3.787  94.189  60.580  1.00 23.29           H   new
ATOM      0  HB3 ARG B   5       2.773  95.204  61.586  1.00 23.29           H   new
ATOM      0  HG2 ARG B   5       4.638  95.959  62.907  1.00 23.24           H   new
ATOM      0  HG3 ARG B   5       5.679  94.674  62.330  1.00 23.24           H   new
ATOM      0  HD2 ARG B   5       5.701  95.615  60.091  1.00 29.00           H   new
ATOM      0  HD3 ARG B   5       4.512  96.821  60.541  1.00 29.00           H   new
ATOM      0  HE  ARG B   5       6.199  98.190  61.259  1.00 32.94           H   new
ATOM      0 HH11 ARG B   5       7.504  94.867  61.294  1.00 30.53           H   new
ATOM      0 HH12 ARG B   5       9.138  95.423  61.671  1.00 30.53           H   new
ATOM      0 HH21 ARG B   5       8.298  98.855  61.593  1.00 35.22           H   new
ATOM      0 HH22 ARG B   5       9.585  97.670  61.837  1.00 35.22           H   new
ATOM    112  N   GLN B   9      10.186  91.359  59.521  1.00 15.79           N
ATOM    113  CA  GLN B   9       9.578  92.637  59.819  1.00 18.37           C
ATOM    114  C   GLN B   9       8.076  92.653  59.539  1.00 17.17           C
ATOM    115  O   GLN B   9       7.278  93.226  60.303  1.00 16.99           O
ATOM    116  CB  GLN B   9      10.229  93.787  59.068  1.00 19.44           C
ATOM    117  CG  GLN B   9      11.591  94.142  59.641  1.00 24.67           C
ATOM    118  CD  GLN B   9      11.455  94.628  61.084  0.00 27.67           C
ATOM    119  NE2 GLN B   9      11.988  94.028  62.025  1.00 31.74           O   flip
ATOM    120  OE1 GLN B   9      10.723  95.704  61.263  1.00 31.80           N   flip
ATOM      0  HA  GLN B   9       9.716  92.803  60.887  1.00 18.37           H   new
ATOM      0  HB2 GLN B   9      10.336  93.519  58.017  1.00 19.44           H   new
ATOM      0  HB3 GLN B   9       9.579  94.661  59.109  1.00 19.44           H   new
ATOM      0  HG2 GLN B   9      12.245  93.271  59.604  1.00 24.67           H   new
ATOM      0  HG3 GLN B   9      12.057  94.917  59.032  1.00 24.67           H   new
ATOM      0 HE21 GLN B   9      12.550  93.196  61.846  1.00 31.74           H   new
ATOM      0 HE22 GLN B   9      11.866  94.365  62.980  1.00 31.74           H   new
HETATM 7573  O   HOH B 407       6.633  88.635  66.643  1.00 11.41           O
HETATM 7693  O   HOH B 527       3.491  85.743  61.660  1.00  4.12           O
'''

#theoretically this should flip when symmetry is detected
sym_flip = '''
CRYST1   14.000   11.000   14.000  90.00  90.00  90.00 P 21 21 21    8
ATOM   1837  N   HIS B   4       0.0    -1.090   1.356  1.00  2.00           N
ATOM   1838  CA  HIS B   4       0.0    -0.537   0.00   1.00  2.00           C
ATOM   1839  C   HIS B   4       0.0     1.000   0.0    1.00  2.00           C
ATOM   1840  O   HIS B   4       0.588   1.645  -0.904  1.00  2.00           O
ATOM   1841  CB  HIS B   4      -1.212  -1.046  -0.789  1.00  2.00           C
ATOM   1842  CG  HIS B   4      -1.164  -2.506  -1.107  1.00  2.00           C
ATOM   1843  ND1 HIS B   4      -2.193  -3.37   -0.7915 1.00  2.00           N
ATOM   1844  CD2 HIS B   4      -0.225  -3.252  -1.734  1.00  2.00           C
ATOM   1845  CE1 HIS B   4      -1.880  -4.586  -1.198  1.00  2.00           C
ATOM   1846  NE2 HIS B   4      -0.696  -4.541  -1.781  1.00  2.00           N
HETATM 7693  O   HOH B 527      -2.263  -7.663  -4.291  1.00  4.12           O
'''

class TestReduce(unittest.TestCase) :

  def test_default_reduce(self) :
    args = [reduce_exe,'-q','-DB',het_dict_path,'-']
    pop = subprocess.Popen(args,stdout=subprocess.PIPE,stdin=subprocess.PIPE)
    out = pop.communicate(input=default_pdb_str)[0]
    #fle = open('default_reduce.pdb','w')
    #fle.write(out)
    #fle.close()
    self.assertEqual(out,default_answer)

  def test_flip_reduce(self) :
    args = [reduce_exe,'-q','-flip','-DB',het_dict_path,'-']
    pop = subprocess.Popen(args,stdout=subprocess.PIPE,stdin=subprocess.PIPE)
    out = pop.communicate(input=flip_test_str)[0]
    #fle = open('flip_reduce.pdb','w')
    #fle.write(out)
    #fle.close()
    self.assertEqual(out,flip_answer)

  def test_nuc_reduce(self) :
    args = [reduce_exe,'-q','-nuc','-flip','-DB',het_dict_path,'-']
    pop = subprocess.Popen(args,stdout=subprocess.PIPE,stdin=subprocess.PIPE)
    out = pop.communicate(input=flip_test_str)[0]
    #fle = open('nuc_reduce.pdb','w')
    #fle.write(out)
    #fle.close()
    self.assertEqual(out,nuc_answer)
    self.assertNotEqual(out,flip_answer)

if __name__ == '__main__' :
    unittest.main()
#suite = unittest.TestLoader().loadTestsFromTestCase(TestReduce)
#unittest.TextTestRunner(verbosity=2).run(suite)
