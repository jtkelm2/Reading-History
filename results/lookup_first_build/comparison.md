| Pattern | search | tables |
|---|---|---|
| 001_traditional_Sailboat | cost 0, 0.0 s | cost 0, 0.1 s |
| 002_traditional_Kabuto | cost 0, 0.0 s | cost 0, 0.0 s |
| 004_traditional_Crane | cost 4, 0.0 s | cost 4, 0.1 s |
| 006_ku_Thirds_Pinwheel | cost 4, 0.0 s | cost 4, 0.1 s |
| 007_ku_Pinwheel_Pockets | cost 0, 0.0 s | cost 0, 0.0 s |
| 015_ku_Iso-Area_Half_Square | cost 610*, 35.2 s | no local history (degree 6) |
| 016_tanakaoripa_enarc-crane_final | cost 13*, 57.8 s | cost 15*, 4.9 s |
| 018_komatsu_Butterfly_Base | timeout | no local history (degree 8, not tabulated) |
| 026_furutaoripa_orihazuru | cost 0, 7.3 s | no local history (degree 16, not tabulated) |
| 030_tanakaoripa_fish_base | timeout | no local history (degree 8, not tabulated) |
| 035_tachioripa_kamehameha | cost 3*, 0.1 s | cost 3*, 0.3 s |
| 038_traditionaloripa_Yakko | cost 6, 0.0 s | cost 6, 0.1 s |
| 040_komatsu_Gentle_Dragon | cost 24, 2.4 s | no local history (degree 8, not tabulated) |
| 043_ku_Quad_Overlap | cost 0, 0.0 s | cost 0, 0.0 s |
| 045_traditionaloripa_Pig | cost 4, 0.1 s | no local history (degree 8, not tabulated) |
| 048_komatsu_Dutch_Rabbit | cost 19*, 25.9 s | cost 25*, 20.5 s |
| 052_traditionaloripa_House | cost 2, 0.0 s | cost 2, 0.1 s |
| 054_tanakaoripa_Hexagon | cost 605*, 18.2 s | no local history (degree 4) |
| 055_komatsu_Sheep | timeout | no local history (degree 12, not tabulated) |
| 056_tachioripa_Hoodman | cost 2, 0.0 s | cost 2, 0.1 s |
| 058_ku_Color-change_Pinwheel_1 | cost 16*, 0.1 s | cost 17*, 1.8 s |
| 062_furutaoripa_houou | timeout | no local history (degree 20, not tabulated) |
| 067_traditionaloripa_damashi_bune | cost 2, 0.0 s | cost 2, 0.1 s |
| 070_ku_Pygmy_Jerboa | cost 9*, 13.3 s | cost 13*, 5.7 s |
| 073_komatsu_Squirrel | cost 22, 1.6 s | cost 22, 0.4 s |
| 074_miuraoripa_Miura-ori | cost 0, 0.0 s | cost 0, 0.1 s |
| 082_traditionaloripa_4_Birdbase | cost 16, 0.9 s | no local history (degree 8, not tabulated) |
| 084_ku_4x4_Checkerboard | cost 614*, 17.8 s | no local history (degree 4) |
| 089_traditionaloripa_9_Birdbase | timeout | no local history (degree 8, not tabulated) |
| 091_ku_4_Stripes_1 | cost 6*, 0.3 s | no local history (degree 8, not tabulated) |
| 093_ku_4_Stripes_3 | cost 5, 0.2 s | no local history (degree 6) |
| 096_komatsu_Horned_Owl | cost 16*, 0.6 s | no local history (degree 8, not tabulated) |
| 102_traditionaloripa_Pinwheel | cost 4, 0.0 s | cost 4, 0.1 s |
| 113_komatsu_Fox | cost 8, 0.1 s | no local history (degree 8, not tabulated) |
| 114_kei_Kirby | cost 16, 0.3 s | cost 18*, 0.7 s |
| 141_kei_kitsune | cost 17, 0.5 s | cost 17, 0.2 s |
| 186_kei_Sus | cost 8*, 0.6 s | no local history (degree 8, not tabulated) |
| 244_komatsuoripa_Smallbird_Base | cost 2, 0.0 s | cost 2, 0.1 s |
| 292_tanaka_Aperiodic_Monotile_1 | cost 5*, 0.0 s | cost 5*, 0.1 s |
| 301_boxhard_Assigned_Crossover | cost 0, 0.0 s | cost 0, 0.0 s |
| 302_boxhard_Assigned_Clause | cost 16*, 7.0 s | cost 16*, 3.9 s |
| 303_boxhard_Assigned_Split | cost 6*, 0.3 s | cost 6*, 0.2 s |
| 304_boxhard_Assigned_Solvable | timeout | timeout |
| 315_kei_Marnie_(Pokemon) | timeout | no local history (degree 6) |
| 316_boxhard_Unassigned_Clause | cost 1, 0.0 s | cost 1, 0.1 s |
| 317_boxhard_Unassigned_Crossover | cost 0, 0.0 s | cost 0, 0.0 s |
| 318_boxhard_Unassigned_Split | cost 2, 0.0 s | cost 2, 0.0 s |
| 325_bernhayes_Unassigned_45˚_Crossover | cost 1200*, 0.6 s | no local history (degree 4) |

resolved: search 41, tables 28, both 28, only search 13, only tables 0
on patterns resolved by both: same cost 23 of 28; search cheaper 5; tables cheaper 0
search time on those: total 106.9 s, median 0.00 s, max 57.8 s
tables time on those: total 39.8 s, median 0.10 s, max 20.5 s
tables failures by reason:
  no local history (degree 12, not tabulated): 1  055_komatsu_Sheep
  no local history (degree 16, not tabulated): 1  026_furutaoripa_orihazuru
  no local history (degree 20, not tabulated): 1  062_furutaoripa_houou
  no local history (degree 4): 3  054_tanakaoripa_Hexagon, 084_ku_4x4_Checkerboard, 325_bernhayes_Unassigned_45˚_Crossover
  no local history (degree 6): 3  015_ku_Iso-Area_Half_Square, 093_ku_4_Stripes_3, 315_kei_Marnie_(Pokemon)
  no local history (degree 8, not tabulated): 10  018_komatsu_Butterfly_Base, 030_tanakaoripa_fish_base, 040_komatsu_Gentle_Dragon, 045_traditionaloripa_Pig, 082_traditionaloripa_4_Birdbase, 089_traditionaloripa_9_Birdbase, 091_ku_4_Stripes_1, 096_komatsu_Horned_Owl, 113_komatsu_Fox, 186_kei_Sus
  timeout: 1  304_boxhard_Assigned_Solvable
