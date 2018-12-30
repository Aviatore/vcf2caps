# VCF2CAPS v2.0 - the software for CAPS markers identification from Variant Call Format (VCF) file.
# Copyright 2018 Wojciech Wesołowski

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>

use strict;
use warnings;

use Tk;
use Tk::NoteBook;
use Tk::FBox;
use Tk::Animation;
use Tk::Listbox;
use Tk::Pane;
use Tk::ProgressBar;
use Tk::PNG;

use threads;
use threads::shared;
use LWP::Simple;
use Digest::MD5 qw(md5_hex);
use Data::Dumper;
use Encode;
use utf8;

$|++;

my $actualSNPNo:shared = 0;
my $caps_filtered:shared = 0;
my $capsMining_percent;
my $cfw_c2f_convertion_percent;
my $cfw_c2f_input_file:shared = "";
my $cfw_c2f_output_file:shared = "";
my $cfw_gf_CAPS_filtering_percent;
my $cfw_gf_group_1_maxError_value = 0;
my $cfw_gf_group_2_maxError_value = 0;
my $cfw_gf_group_3_maxError_value = 0;
my $cfw_gf_input_file:shared = "";
my $cfw_gf_output_file:shared = "";
my $cfw_scf_CAPS_filtering_percent;
my $cfw_scf_input_file:shared = "";
my $cfw_scf_output_file:shared = "";
my $comp_state = 0;
my $custom = 0;
my $custom_value = "";
my $cutters4 = 0;
my $cutters5 = 0;
my $cutters6 = 0;
my $die:shared = 0;
my $die_confirm:shared = 0;
my $download_enzyme_db_result:shared = 0;
my $enzyme_check;
my $enzyme_check_status;
my $enzyme_file_name:shared;
my $enzyme_start_analysis:shared = 0;
my $enzyme_zip:shared = 0;
my $enzyme_zip_stop:shared = 0;
my $fh;
my $genome_exists = 0;
my $iso_state = 0;
my $jobID:shared = 0;
my $L_center_col2_1_entry;
my $L_center_col2_2_entry;
my $line_vcf:shared = 0;
my $log_first_use = 0;
my $numberOfSNPs:shared = 0;
my $numberOfSNPsAfter:shared;
my $numberOfSNPsBefore:shared = 0;
my $oneFiltGroup_twoGenotypes:shared;
my $onlyCommercially = 1;
my $out;
my $output_seq_len:shared = 500;
my $raw_vcf_check;
my $raw_vcf_check_status;
my $raw_vcf_file_name:shared;
my $reference_check;
my $reference_check_status;
my $reference_file_name:shared;
my $reference_md5:shared = "";
my $reference_start_analysis:shared = 0;
my $selected_enz_name;
my $singleCutSite_percent;
my $snps_seq_len:shared = 40;
my $snpsNo:shared = 0;
my $stop:shared;
my $sVCF_file_name:shared;
my $terminal;
my $total_caps_number:shared = 0;
my $v2c_file_name:shared = "";
my $vcf_md5:shared = "";
my $vcf_start_analysis:shared = 0;
my $working_dir:shared = "";
my %cfw_groups:shared;
$cfw_groups{1} = &share({});
$cfw_groups{2} = &share({});
$cfw_groups{3} = &share({});
my %enzymes_db:shared;
my %enzymes_for_analysis:shared;
my %raw_vcf_analysis_results:shared = (err_code => 0);
my %reference_index_check_result:shared = ();
my %v2c_check_result:shared = (err_type => "");
my %vcf_analysis_results:shared = (err_code => 0);
my %vcf2capsOutput_results:shared = (err_code => 0);
my @allEnzymesNames:shared;
my @caps_filtration_result:shared = ();
my @caps_mining_results:shared;
my @caps_to_fasta_result:shared = ();
my @enzyme_analysis_results:shared = (0);
my @linie:shared;
my @markersOnTheEdge:shared = ();
my @reference_analysis_results:shared = (0);
my @selected_enz_names:shared = ();
my @seqExtractor_error_code:shared = ("",0);
my @sequencesNotPresentInRef:shared = ();
my @singleCutSite_results:shared;


print 'VCF2CAPS v2.0 - the software for CAPS markers identification 
from Variant Call Format (VCF) file.
Copyright 2018 Wojciech Wesolowski

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>

This window is required for VCF2CAPS to run properly.
You can minimise it but please do not close it.

';

#------------------------------------------------------------------------------#
# Detaching the subroutine 'work', that will process all analysis.             #
# The 'work' subroutine runs in the background as a loop, checking the value   #
# of the $jobID variable that determines which analysis should be performed.   #
# Thanks to assigning computational demanding processes to the new subroutine  #
# running in a new thread, the main window does not freeze and other tasks     #
# can be performed in the meantime. Perl Tk module is not thread safe and any  #
# attempt to create a new thread within Tk raise an error. The below (probably #
# odd) solution was my answer to that problem.                                 #
#------------------------------------------------------------------------------#
threads->create( \&work )->detach();

my $mw = MainWindow->new();
$mw->withdraw;
$mw->title("VCF2CAPS v2.0");
$mw->minsize( qw(500 300) );


#-----------------------#
# Images base64 encoded #
#-----------------------#
my $warning_image = $mw->Photo(-data => 'R0lGODlhEAAQAOZRAP/eAPPFc/alEu7s6Py+AfarJfzZm/LPkP7shv/iAP/mAO/iy/aqIPm4OPzV
fP/4mfanF//jKvejC/7vtv/wY//zlv/3dP/4d/73e/vPdf/vWf/wq//kLfzSc//tKf/4Yvq9Gv/j
B//plf/bAv/cMv/3hvHZrfinFvTFcfq8GPvCRPfBXf3prPPHePm2M/vQeP3lnvPGdvHZr/fBXvzQ
VP/nivzUZP/qj//YBvTGcvS7VP/tNP/cNP7urf7rpv/kEf7obP/ZC/rAQP/3tv7nY//5rfW7Vfin
F//4XP3mif/XAP/TAPLKgf/bAP/smvegA////+7u7gAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACH/
C1hNUCBEYXRhWE1QPD94cGFja2V0IGJlZ2luPSLvu78iIGlkPSJXNU0wTXBDZWhpSHpyZVN6TlRj
emtjOWQiPz4gPHg6eG1wbWV0YSB4bWxuczp4PSJhZG9iZTpuczptZXRhLyIgeDp4bXB0az0iQWRv
YmUgWE1QIENvcmUgNS4wLWMwNjAgNjEuMTM0Nzc3LCAyMDEwLzAyLzEyLTE3OjMyOjAwICAgICAg
ICAiPiA8cmRmOlJERiB4bWxuczpyZGY9Imh0dHA6Ly93d3cudzMub3JnLzE5OTkvMDIvMjItcmRm
LXN5bnRheC1ucyMiPiA8cmRmOkRlc2NyaXB0aW9uIHJkZjphYm91dD0iIiB4bWxuczp4bXA9Imh0
dHA6Ly9ucy5hZG9iZS5jb20veGFwLzEuMC8iIHhtbG5zOnhtcE1NPSJodHRwOi8vbnMuYWRvYmUu
Y29tL3hhcC8xLjAvbW0vIiB4bWxuczpzdFJlZj0iaHR0cDovL25zLmFkb2JlLmNvbS94YXAvMS4w
L3NUeXBlL1Jlc291cmNlUmVmIyIgeG1wOkNyZWF0b3JUb29sPSJBZG9iZSBQaG90b3Nob3AgQ1M1
IFdpbmRvd3MiIHhtcE1NOkluc3RhbmNlSUQ9InhtcC5paWQ6OEYyMkVFOTFFMzcxMTFFODgwNzJE
MDA0QTMyMUExMUEiIHhtcE1NOkRvY3VtZW50SUQ9InhtcC5kaWQ6OEYyMkVFOTJFMzcxMTFFODgw
NzJEMDA0QTMyMUExMUEiPiA8eG1wTU06RGVyaXZlZEZyb20gc3RSZWY6aW5zdGFuY2VJRD0ieG1w
LmlpZDo4RjIyRUU4RkUzNzExMUU4ODA3MkQwMDRBMzIxQTExQSIgc3RSZWY6ZG9jdW1lbnRJRD0i
eG1wLmRpZDo4RjIyRUU5MEUzNzExMUU4ODA3MkQwMDRBMzIxQTExQSIvPiA8L3JkZjpEZXNjcmlw
dGlvbj4gPC9yZGY6UkRGPiA8L3g6eG1wbWV0YT4gPD94cGFja2V0IGVuZD0iciI/PgH//v38+/r5
+Pf29fTz8vHw7+7t7Ovq6ejn5uXk4+Lh4N/e3dzb2tnY19bV1NPS0dDPzs3My8rJyMfGxcTDwsHA
v769vLu6ubi3trW0s7KxsK+urayrqqmop6alpKOioaCfnp2cm5qZmJeWlZSTkpGQj46NjIuKiYiH
hoWEg4KBgH9+fXx7enl4d3Z1dHNycXBvbm1sa2ppaGdmZWRjYmFgX15dXFtaWVhXVlVUU1JRUE9O
TUxLSklIR0ZFRENCQUA/Pj08Ozo5ODc2NTQzMjEwLy4tLCsqKSgnJiUkIyIhIB8eHRwbGhkYFxYV
FBMSERAPDg0MCwoJCAcGBQQDAgEAACH5BAEAAFEALAAAAAAQABAAAAeegFGCg4SFhlFMBQVMh4QB
IEUPKQGNghAYSB8lApVMNhdPEhY0jIYLDRVHUFAnQy4Lhig9Hk+rTztJOYUHDhoKtVBPChQdB4MD
KwghCcBPCT8IMwOCLTAcAADN2BEsMVEmL0RN400EBOQjQBkyRhNBSvDN8Eo4PjoCGyRL+837Szw3
BARQUcOJQScGDBx0IkIIJSYMnkicSJFBqUqGAgEAOw==');

my $icon_base64 = 'iVBORw0KGgoAAAANSUhEUgAAACAAAAAhCAYAAAC4JqlRAAAACXBIWXMAAAsTAAALEwEAmpwYAAAK
T2lDQ1BQaG90b3Nob3AgSUNDIHByb2ZpbGUAAHjanVNnVFPpFj333vRCS4iAlEtvUhUIIFJCi4AU
kSYqIQkQSoghodkVUcERRUUEG8igiAOOjoCMFVEsDIoK2AfkIaKOg6OIisr74Xuja9a89+bN/rXX
Pues852zzwfACAyWSDNRNYAMqUIeEeCDx8TG4eQuQIEKJHAAEAizZCFz/SMBAPh+PDwrIsAHvgAB
eNMLCADATZvAMByH/w/qQplcAYCEAcB0kThLCIAUAEB6jkKmAEBGAYCdmCZTAKAEAGDLY2LjAFAt
AGAnf+bTAICd+Jl7AQBblCEVAaCRACATZYhEAGg7AKzPVopFAFgwABRmS8Q5ANgtADBJV2ZIALC3
AMDOEAuyAAgMADBRiIUpAAR7AGDIIyN4AISZABRG8lc88SuuEOcqAAB4mbI8uSQ5RYFbCC1xB1dX
Lh4ozkkXKxQ2YQJhmkAuwnmZGTKBNA/g88wAAKCRFRHgg/P9eM4Ors7ONo62Dl8t6r8G/yJiYuP+
5c+rcEAAAOF0ftH+LC+zGoA7BoBt/qIl7gRoXgugdfeLZrIPQLUAoOnaV/Nw+H48PEWhkLnZ2eXk
5NhKxEJbYcpXff5nwl/AV/1s+X48/Pf14L7iJIEyXYFHBPjgwsz0TKUcz5IJhGLc5o9H/LcL//wd
0yLESWK5WCoU41EScY5EmozzMqUiiUKSKcUl0v9k4t8s+wM+3zUAsGo+AXuRLahdYwP2SycQWHTA
4vcAAPK7b8HUKAgDgGiD4c93/+8//UegJQCAZkmScQAAXkQkLlTKsz/HCAAARKCBKrBBG/TBGCzA
BhzBBdzBC/xgNoRCJMTCQhBCCmSAHHJgKayCQiiGzbAdKmAv1EAdNMBRaIaTcA4uwlW4Dj1wD/ph
CJ7BKLyBCQRByAgTYSHaiAFiilgjjggXmYX4IcFIBBKLJCDJiBRRIkuRNUgxUopUIFVIHfI9cgI5
h1xGupE7yAAygvyGvEcxlIGyUT3UDLVDuag3GoRGogvQZHQxmo8WoJvQcrQaPYw2oefQq2gP2o8+
Q8cwwOgYBzPEbDAuxsNCsTgsCZNjy7EirAyrxhqwVqwDu4n1Y8+xdwQSgUXACTYEd0IgYR5BSFhM
WE7YSKggHCQ0EdoJNwkDhFHCJyKTqEu0JroR+cQYYjIxh1hILCPWEo8TLxB7iEPENyQSiUMyJ7mQ
AkmxpFTSEtJG0m5SI+ksqZs0SBojk8naZGuyBzmULCAryIXkneTD5DPkG+Qh8lsKnWJAcaT4U+Io
UspqShnlEOU05QZlmDJBVaOaUt2ooVQRNY9aQq2htlKvUYeoEzR1mjnNgxZJS6WtopXTGmgXaPdp
r+h0uhHdlR5Ol9BX0svpR+iX6AP0dwwNhhWDx4hnKBmbGAcYZxl3GK+YTKYZ04sZx1QwNzHrmOeZ
D5lvVVgqtip8FZHKCpVKlSaVGyovVKmqpqreqgtV81XLVI+pXlN9rkZVM1PjqQnUlqtVqp1Q61Mb
U2epO6iHqmeob1Q/pH5Z/YkGWcNMw09DpFGgsV/jvMYgC2MZs3gsIWsNq4Z1gTXEJrHN2Xx2KruY
/R27iz2qqaE5QzNKM1ezUvOUZj8H45hx+Jx0TgnnKKeX836K3hTvKeIpG6Y0TLkxZVxrqpaXllir
SKtRq0frvTau7aedpr1Fu1n7gQ5Bx0onXCdHZ4/OBZ3nU9lT3acKpxZNPTr1ri6qa6UbobtEd79u
p+6Ynr5egJ5Mb6feeb3n+hx9L/1U/W36p/VHDFgGswwkBtsMzhg8xTVxbzwdL8fb8VFDXcNAQ6Vh
lWGX4YSRudE8o9VGjUYPjGnGXOMk423GbcajJgYmISZLTepN7ppSTbmmKaY7TDtMx83MzaLN1pk1
mz0x1zLnm+eb15vft2BaeFostqi2uGVJsuRaplnutrxuhVo5WaVYVVpds0atna0l1rutu6cRp7lO
k06rntZnw7Dxtsm2qbcZsOXYBtuutm22fWFnYhdnt8Wuw+6TvZN9un2N/T0HDYfZDqsdWh1+c7Ry
FDpWOt6azpzuP33F9JbpL2dYzxDP2DPjthPLKcRpnVOb00dnF2e5c4PziIuJS4LLLpc+Lpsbxt3I
veRKdPVxXeF60vWdm7Obwu2o26/uNu5p7ofcn8w0nymeWTNz0MPIQ+BR5dE/C5+VMGvfrH5PQ0+B
Z7XnIy9jL5FXrdewt6V3qvdh7xc+9j5yn+M+4zw33jLeWV/MN8C3yLfLT8Nvnl+F30N/I/9k/3r/
0QCngCUBZwOJgUGBWwL7+Hp8Ib+OPzrbZfay2e1BjKC5QRVBj4KtguXBrSFoyOyQrSH355jOkc5p
DoVQfujW0Adh5mGLw34MJ4WHhVeGP45wiFga0TGXNXfR3ENz30T6RJZE3ptnMU85ry1KNSo+qi5q
PNo3ujS6P8YuZlnM1VidWElsSxw5LiquNm5svt/87fOH4p3iC+N7F5gvyF1weaHOwvSFpxapLhIs
OpZATIhOOJTwQRAqqBaMJfITdyWOCnnCHcJnIi/RNtGI2ENcKh5O8kgqTXqS7JG8NXkkxTOlLOW5
hCepkLxMDUzdmzqeFpp2IG0yPTq9MYOSkZBxQqohTZO2Z+pn5mZ2y6xlhbL+xW6Lty8elQfJa7OQ
rAVZLQq2QqboVFoo1yoHsmdlV2a/zYnKOZarnivN7cyzytuQN5zvn//tEsIS4ZK2pYZLVy0dWOa9
rGo5sjxxedsK4xUFK4ZWBqw8uIq2Km3VT6vtV5eufr0mek1rgV7ByoLBtQFr6wtVCuWFfevc1+1d
T1gvWd+1YfqGnRs+FYmKrhTbF5cVf9go3HjlG4dvyr+Z3JS0qavEuWTPZtJm6ebeLZ5bDpaql+aX
Dm4N2dq0Dd9WtO319kXbL5fNKNu7g7ZDuaO/PLi8ZafJzs07P1SkVPRU+lQ27tLdtWHX+G7R7ht7
vPY07NXbW7z3/T7JvttVAVVN1WbVZftJ+7P3P66Jqun4lvttXa1ObXHtxwPSA/0HIw6217nU1R3S
PVRSj9Yr60cOxx++/p3vdy0NNg1VjZzG4iNwRHnk6fcJ3/ceDTradox7rOEH0x92HWcdL2pCmvKa
RptTmvtbYlu6T8w+0dbq3nr8R9sfD5w0PFl5SvNUyWna6YLTk2fyz4ydlZ19fi753GDborZ752PO
32oPb++6EHTh0kX/i+c7vDvOXPK4dPKy2+UTV7hXmq86X23qdOo8/pPTT8e7nLuarrlca7nuer21
e2b36RueN87d9L158Rb/1tWeOT3dvfN6b/fF9/XfFt1+cif9zsu72Xcn7q28T7xf9EDtQdlD3YfV
P1v+3Njv3H9qwHeg89HcR/cGhYPP/pH1jw9DBY+Zj8uGDYbrnjg+OTniP3L96fynQ89kzyaeF/6i
/suuFxYvfvjV69fO0ZjRoZfyl5O/bXyl/erA6xmv28bCxh6+yXgzMV70VvvtwXfcdx3vo98PT+R8
IH8o/2j5sfVT0Kf7kxmTk/8EA5jz/GMzLdsAAAAgY0hSTQAAeiUAAICDAAD5/wAAgOkAAHUwAADq
YAAAOpgAABdvkl/FRgAABv9JREFUeNrEl21MlNkVx3/PMzMOIFRGQSgvGqMYXqRKLKsGUIxaY0gQ
azURtW6oWSVxm64mpmmbfukHywebdY1tUIdk0UWjYUlMSwxG5C1EcN1AIIrCyBAsLzO8TBw6DDPz
PLcfWK4Ob6vJJp7kSe5z7jnn/u+555x7riKEEHxEUvfs2YPZbGb58uXyCw8PJy0tDbvd/pMsUlxc
TFhYWNAaZrOZ3bt3g9frFUeOHBGA/FasWCF6enrET0kHDx4MWqOwsFB4vV6BEEL09/eL6OhoORkX
FydcLtePGtU0TbjdbjEyMiKGh4eFw+EQY2NjYnJyco7siRMngjbY19cnhBDCCJCQkEBRURElJSUA
DAwMUFNTw6FDh+a40+Vy0d3djd1uZ3h4GI/Hg8/nQ9M0FEXBaDRiNpuJiIggJiaGlJQUjEYjNTU1
0sbx48dZtWoVAMpMEPb09JCRkcHExAQA+/bto7q6WioNDw/T1NREV1cXfr9/WllRfvT8FUVBVVU6
Ozt58uQJTqeTzs5O4uPjURTlLYAZZDdv3gTAaDTS2tpKRkYGDQ0NNDQ0oGnajFXMoUuJj/05MbEr
iYiIYMkSM0LX8Xo9uFwunE4nIyMj6LqOrusoikIgEGB8fJwDBw7gcDjIy8sLBtDU1EROTo5Ef+bM
GXbu3ElHR8d0yqgqP1seRdYnvyQ1NQVVASEEIWHh83pjcnISq9XK06dPSUpKQgiBoijouk5eXh4Z
GRkY31XIzs4mNzeXuro6AK5fv47ZbCYyMhIBJK1bz68P7MfvneQ///wr7lEnKiqTbje7Pitm/cY9
QQBCQ0Opra2lqqqK5ORk9u7dS1RUFImJiWzatGlaaHa0Wq3WoHTZtWuXKCkpEa2t3wkhhHD8t1tc
/PRXouXf3whNCwghhOhurxMXP90t3owMBdlqaWkRBoNB2rp8+bK4dOmS6O3tlTJzALjdbpGSkiKV
LBaLqKmpkfM+n1e8aH8sdF0L0vvXF78RL9q/C+KdPHlS2tm4caMYGBgQDx8+DJJhvvw+e/ZskBes
Vuui9aDxrlVcO39U+H0+ybPb7cJisUgbFy9enFdXnR04Q0NDmEwmwsPDJe/KlStvM+Adevl9HV//
+XcM93Vx5E9fYTSZ5Fx5eTnj4+MAxMbGUlhYOH+ezkbU3NwsLly4ILZt2xbkhYaGBikz5f2fuPP3
8+Lrv3wm7M+/n7OrN2/eiLVr10rdc+fOLei9OR4YHBxE0zQ2b96M0fg2Sa5duybHlV+eIypxDb/9
WymrkzNmNoL+g5cqKyux2WwyE4qKihYsVMbZDLfbTSAQIC4ujpSUFFkDqqqqGBocwDc2RFdjA6ad
EVR++UeEpoGi4J14w/rNOXySdzQI7P79+0lNTX1/ADNl1mAwsH37dmw2Gx6Ph4mJCa5eu84fzhRz
6h+30HQdoetSTw/4CV8RQ1NjA83NzZJ/+vTpRUv1HAAGgwEATdPIzc2lt7dX3gllZWX8/vPPiVv/
iwUNfvXFeTnOzc1lx44dizcksxkRERHSE6GhoZw6dUrO9fX1Ufnttwsae/bsGffu3Xvv3c8LIDo6
WgaVw+EgJyeHpKSkoGDU33H9u2S1WpmamgIgLS2N/Pz8DwewevVqOX79+jWjo6NBO2lpaeHBgwdz
DDmdTioqKuR/UVERoaGhHw4gMTFRekHXdZ4/f86xY8ewWCxS5urVq3MMVVRUMDQ0tGjh0XWd+vp6
RkZGFgZgMBjIzs5m5pbu6OhgcHCQ4uJiKVNdXS3TE2Bqagqr1Sr/jx49Smxs7NzdqioDAwPU1tYu
DGDm/BISEmQwNjY2smXLFlmevV4vZWVl8wJarPDYbDZ6e3vp6uqira1tYQAGg4H8/HxCQkIAGBsb
o7u7m8zMTClz+/ZtnE4nAKWlpZJfUFAwp/AIIXj8+DG3bt3C7/cjhOD+/fu0trYGd0Szqaenh7t3
7+Lz+VBVFVVVaW9vp66uDqfTSWlpKdnZ2WzYsEEeWVNTE1lZWXLhV69eUV9fT39/v7S7Zs0a0tPT
GR0dXRzATCZUVlbicrkAMJvNTE5OYrPZWLJkCUuXLuXGjRv4/X6ysrIoLy/H4/Fgt9t58eIFDodD
gjGZTGRlZZGTk4PBYJhu0d7naebxeHj06BFtbW3ouo7RaERVp0/P5/PJtjwsLAyj0YjX6w3qERVF
ITU1la1bt8rYknMf8jYcGxujpaWF6upqLBYLISEhGAyG6fb6h2ZT0zT8fj+qqhIdHc26detIT08n
JiZm/rb9Qx+nXV1dZGZmYjabiYqKYtmyZYSEhGAymTh8+DDJycnEx8djsViIjIyUnnrvy2gxevny
JQUFBdN3v67jcDjkGQcCAQKBAHfu3GHlypXvbVP52M9zI6B8TAD/HwC4utyA+72HQQAAAABJRU5E
rkJggg==';

my $logo = $mw->Photo(-data => $icon_base64, -format => 'PNG', -width => 32, -height => 31);
$mw->iconimage($logo);

my $folder_image = $mw->Photo(-data => 'R0lGODlhEAAQANUyAPT09Obm5vDw8MC6r+fi18nJyezs7PX19a+vr7awperq6t3XzETXAPLy8ujo
6J7lAIjhAM/Kv93Vxu3t7ba2trW1tbToAPf39xfQALu7u8rKyrq6uufh17i4uNrSwy3TANbW1vj4
+MO+s+vr69TPxMnDuNnTyNvUxfn5+eDg4Nvb28XFxXHeAFraAL29vaCgoPr6+v///+7u7gAAAAAA
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACH/C1hNUCBEYXRhWE1QPD94cGFja2V0
IGJlZ2luPSLvu78iIGlkPSJXNU0wTXBDZWhpSHpyZVN6TlRjemtjOWQiPz4gPHg6eG1wbWV0YSB4
bWxuczp4PSJhZG9iZTpuczptZXRhLyIgeDp4bXB0az0iQWRvYmUgWE1QIENvcmUgNS4wLWMwNjAg
NjEuMTM0Nzc3LCAyMDEwLzAyLzEyLTE3OjMyOjAwICAgICAgICAiPiA8cmRmOlJERiB4bWxuczpy
ZGY9Imh0dHA6Ly93d3cudzMub3JnLzE5OTkvMDIvMjItcmRmLXN5bnRheC1ucyMiPiA8cmRmOkRl
c2NyaXB0aW9uIHJkZjphYm91dD0iIiB4bWxuczp4bXA9Imh0dHA6Ly9ucy5hZG9iZS5jb20veGFw
LzEuMC8iIHhtbG5zOnhtcE1NPSJodHRwOi8vbnMuYWRvYmUuY29tL3hhcC8xLjAvbW0vIiB4bWxu
czpzdFJlZj0iaHR0cDovL25zLmFkb2JlLmNvbS94YXAvMS4wL3NUeXBlL1Jlc291cmNlUmVmIyIg
eG1wOkNyZWF0b3JUb29sPSJBZG9iZSBQaG90b3Nob3AgQ1M1IFdpbmRvd3MiIHhtcE1NOkluc3Rh
bmNlSUQ9InhtcC5paWQ6NDYwNDUzMDMxMjVCMTFFODg5NzY5MjU5MkU5NEZBODIiIHhtcE1NOkRv
Y3VtZW50SUQ9InhtcC5kaWQ6NDYwNDUzMDQxMjVCMTFFODg5NzY5MjU5MkU5NEZBODIiPiA8eG1w
TU06RGVyaXZlZEZyb20gc3RSZWY6aW5zdGFuY2VJRD0ieG1wLmlpZDo0NjA0NTMwMTEyNUIxMUU4
ODk3NjkyNTkyRTk0RkE4MiIgc3RSZWY6ZG9jdW1lbnRJRD0ieG1wLmRpZDo0NjA0NTMwMjEyNUIx
MUU4ODk3NjkyNTkyRTk0RkE4MiIvPiA8L3JkZjpEZXNjcmlwdGlvbj4gPC9yZGY6UkRGPiA8L3g6
eG1wbWV0YT4gPD94cGFja2V0IGVuZD0iciI/PgH//v38+/r5+Pf29fTz8vHw7+7t7Ovq6ejn5uXk
4+Lh4N/e3dzb2tnY19bV1NPS0dDPzs3My8rJyMfGxcTDwsHAv769vLu6ubi3trW0s7KxsK+urayr
qqmop6alpKOioaCfnp2cm5qZmJeWlZSTkpGQj46NjIuKiYiHhoWEg4KBgH9+fXx7enl4d3Z1dHNy
cXBvbm1sa2ppaGdmZWRjYmFgX15dXFtaWVhXVlVUU1JRUE9OTUxLSklIR0ZFRENCQUA/Pj08Ozo5
ODc2NTQzMjEwLy4tLCsqKSgnJiUkIyIhIB8eHRwbGhkYFxYVFBMSERAPDg0MCwoJCAcGBQQDAgEA
ACH5BAEAADIALAAAAAAQABAAAAaXQJkM4SoWM52KZmIQCl2waBSgiFFWE6cMKoU1AjGrRsuVAhzh
WIXchV0EjoCiw+6iLoAGIMMm+DkcaTEuWi8hEhInHh4xFhYxL4UHkwIGBjEPD5CFAgsLghAQgkIv
IyYRMSyqq6oxpAEkJTEttLW0rjIvKREiggwMo7kqA8QJCTEfH5tOLyAFzxsbMRgYy0IIL9nagghC
QQA7');

my $analyze_image = $mw->Photo(-data => 'R0lGODlhEAAQALMPAACd0QCq4wCv6UeIrwCZzACazwCb0gCh1gCh2ACYzQCk2wCUxgCd0gCe1QBd
lO7u7iH/C1hNUCBEYXRhWE1QPD94cGFja2V0IGJlZ2luPSLvu78iIGlkPSJXNU0wTXBDZWhpSHpy
ZVN6TlRjemtjOWQiPz4gPHg6eG1wbWV0YSB4bWxuczp4PSJhZG9iZTpuczptZXRhLyIgeDp4bXB0
az0iQWRvYmUgWE1QIENvcmUgNS4wLWMwNjAgNjEuMTM0Nzc3LCAyMDEwLzAyLzEyLTE3OjMyOjAw
ICAgICAgICAiPiA8cmRmOlJERiB4bWxuczpyZGY9Imh0dHA6Ly93d3cudzMub3JnLzE5OTkvMDIv
MjItcmRmLXN5bnRheC1ucyMiPiA8cmRmOkRlc2NyaXB0aW9uIHJkZjphYm91dD0iIiB4bWxuczp4
bXA9Imh0dHA6Ly9ucy5hZG9iZS5jb20veGFwLzEuMC8iIHhtbG5zOnhtcE1NPSJodHRwOi8vbnMu
YWRvYmUuY29tL3hhcC8xLjAvbW0vIiB4bWxuczpzdFJlZj0iaHR0cDovL25zLmFkb2JlLmNvbS94
YXAvMS4wL3NUeXBlL1Jlc291cmNlUmVmIyIgeG1wOkNyZWF0b3JUb29sPSJBZG9iZSBQaG90b3No
b3AgQ1M1IFdpbmRvd3MiIHhtcE1NOkluc3RhbmNlSUQ9InhtcC5paWQ6QUEwRUIwNDlEQzZBMTFF
OEJDNzBENTI4MkQ0MTY5RTYiIHhtcE1NOkRvY3VtZW50SUQ9InhtcC5kaWQ6QUEwRUIwNEFEQzZB
MTFFOEJDNzBENTI4MkQ0MTY5RTYiPiA8eG1wTU06RGVyaXZlZEZyb20gc3RSZWY6aW5zdGFuY2VJ
RD0ieG1wLmlpZDpBQTBFQjA0N0RDNkExMUU4QkM3MEQ1MjgyRDQxNjlFNiIgc3RSZWY6ZG9jdW1l
bnRJRD0ieG1wLmRpZDpBQTBFQjA0OERDNkExMUU4QkM3MEQ1MjgyRDQxNjlFNiIvPiA8L3JkZjpE
ZXNjcmlwdGlvbj4gPC9yZGY6UkRGPiA8L3g6eG1wbWV0YT4gPD94cGFja2V0IGVuZD0iciI/PgH/
/v38+/r5+Pf29fTz8vHw7+7t7Ovq6ejn5uXk4+Lh4N/e3dzb2tnY19bV1NPS0dDPzs3My8rJyMfG
xcTDwsHAv769vLu6ubi3trW0s7KxsK+urayrqqmop6alpKOioaCfnp2cm5qZmJeWlZSTkpGQj46N
jIuKiYiHhoWEg4KBgH9+fXx7enl4d3Z1dHNycXBvbm1sa2ppaGdmZWRjYmFgX15dXFtaWVhXVlVU
U1JRUE9OTUxLSklIR0ZFRENCQUA/Pj08Ozo5ODc2NTQzMjEwLy4tLCsqKSgnJiUkIyIhIB8eHRwb
GhkYFxYVFBMSERAPDg0MCwoJCAcGBQQDAgEAACH5BAEAAA8ALAAAAAAQABAAAAQu8MlJq7044+E0
dULnPWA4OkEqZo7RIMp6OUnBHLLlAHyuEwSfbiG0cEbIpHISAQA7');

my $ok_image = $mw->Photo(-data => 'R0lGODlhEAAQALMKAEfbAI7rAC/WAKbwAL31ANX6ABjRAF/hAHbmABizAe7u7gAAAAAAAAAAAAAA
AAAAACH/C1hNUCBEYXRhWE1QPD94cGFja2V0IGJlZ2luPSLvu78iIGlkPSJXNU0wTXBDZWhpSHpy
ZVN6TlRjemtjOWQiPz4gPHg6eG1wbWV0YSB4bWxuczp4PSJhZG9iZTpuczptZXRhLyIgeDp4bXB0
az0iQWRvYmUgWE1QIENvcmUgNS4wLWMwNjAgNjEuMTM0Nzc3LCAyMDEwLzAyLzEyLTE3OjMyOjAw
ICAgICAgICAiPiA8cmRmOlJERiB4bWxuczpyZGY9Imh0dHA6Ly93d3cudzMub3JnLzE5OTkvMDIv
MjItcmRmLXN5bnRheC1ucyMiPiA8cmRmOkRlc2NyaXB0aW9uIHJkZjphYm91dD0iIiB4bWxuczp4
bXA9Imh0dHA6Ly9ucy5hZG9iZS5jb20veGFwLzEuMC8iIHhtbG5zOnhtcE1NPSJodHRwOi8vbnMu
YWRvYmUuY29tL3hhcC8xLjAvbW0vIiB4bWxuczpzdFJlZj0iaHR0cDovL25zLmFkb2JlLmNvbS94
YXAvMS4wL3NUeXBlL1Jlc291cmNlUmVmIyIgeG1wOkNyZWF0b3JUb29sPSJBZG9iZSBQaG90b3No
b3AgQ1M1IFdpbmRvd3MiIHhtcE1NOkluc3RhbmNlSUQ9InhtcC5paWQ6QzE5RDkwMzYxMjVBMTFF
ODgzMEJGOUNFMzY4RjVEOEMiIHhtcE1NOkRvY3VtZW50SUQ9InhtcC5kaWQ6QzE5RDkwMzcxMjVB
MTFFODgzMEJGOUNFMzY4RjVEOEMiPiA8eG1wTU06RGVyaXZlZEZyb20gc3RSZWY6aW5zdGFuY2VJ
RD0ieG1wLmlpZDpDMTlEOTAzNDEyNUExMUU4ODMwQkY5Q0UzNjhGNUQ4QyIgc3RSZWY6ZG9jdW1l
bnRJRD0ieG1wLmRpZDpDMTlEOTAzNTEyNUExMUU4ODMwQkY5Q0UzNjhGNUQ4QyIvPiA8L3JkZjpE
ZXNjcmlwdGlvbj4gPC9yZGY6UkRGPiA8L3g6eG1wbWV0YT4gPD94cGFja2V0IGVuZD0iciI/PgH/
/v38+/r5+Pf29fTz8vHw7+7t7Ovq6ejn5uXk4+Lh4N/e3dzb2tnY19bV1NPS0dDPzs3My8rJyMfG
xcTDwsHAv769vLu6ubi3trW0s7KxsK+urayrqqmop6alpKOioaCfnp2cm5qZmJeWlZSTkpGQj46N
jIuKiYiHhoWEg4KBgH9+fXx7enl4d3Z1dHNycXBvbm1sa2ppaGdmZWRjYmFgX15dXFtaWVhXVlVU
U1JRUE9OTUxLSklIR0ZFRENCQUA/Pj08Ozo5ODc2NTQzMjEwLy4tLCsqKSgnJiUkIyIhIB8eHRwb
GhkYFxYVFBMSERAPDg0MCwoJCAcGBQQDAgEAACH5BAEAAAoALAAAAAAQABAAAAQ0UMlJq70465r2
TEW3JUQpSic6rGkSnG7wUglid/VtJUfP+5gEYDhM7QRI4yVhUAY90KgkAgA7');

my $fail_image = $mw->Photo(-data => 'R0lGODlhEAAQAMQXALgAAMkAAKkAANsAAMAAALAAAJ8AANIAAP8AAKYAANMAAPkAAMYAANkAAKwA
AL8AAN8AAPMAALkAALMAAOwAAOYAAMwAAO7u7gAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACH/C1hN
UCBEYXRhWE1QPD94cGFja2V0IGJlZ2luPSLvu78iIGlkPSJXNU0wTXBDZWhpSHpyZVN6TlRjemtj
OWQiPz4gPHg6eG1wbWV0YSB4bWxuczp4PSJhZG9iZTpuczptZXRhLyIgeDp4bXB0az0iQWRvYmUg
WE1QIENvcmUgNS4wLWMwNjAgNjEuMTM0Nzc3LCAyMDEwLzAyLzEyLTE3OjMyOjAwICAgICAgICAi
PiA8cmRmOlJERiB4bWxuczpyZGY9Imh0dHA6Ly93d3cudzMub3JnLzE5OTkvMDIvMjItcmRmLXN5
bnRheC1ucyMiPiA8cmRmOkRlc2NyaXB0aW9uIHJkZjphYm91dD0iIiB4bWxuczp4bXA9Imh0dHA6
Ly9ucy5hZG9iZS5jb20veGFwLzEuMC8iIHhtbG5zOnhtcE1NPSJodHRwOi8vbnMuYWRvYmUuY29t
L3hhcC8xLjAvbW0vIiB4bWxuczpzdFJlZj0iaHR0cDovL25zLmFkb2JlLmNvbS94YXAvMS4wL3NU
eXBlL1Jlc291cmNlUmVmIyIgeG1wOkNyZWF0b3JUb29sPSJBZG9iZSBQaG90b3Nob3AgQ1M1IFdp
bmRvd3MiIHhtcE1NOkluc3RhbmNlSUQ9InhtcC5paWQ6MjAwRDcyMTQxMjVCMTFFOEE3OUFDRTgy
QTJCN0QxNUIiIHhtcE1NOkRvY3VtZW50SUQ9InhtcC5kaWQ6MjAwRDcyMTUxMjVCMTFFOEE3OUFD
RTgyQTJCN0QxNUIiPiA8eG1wTU06RGVyaXZlZEZyb20gc3RSZWY6aW5zdGFuY2VJRD0ieG1wLmlp
ZDoyMDBENzIxMjEyNUIxMUU4QTc5QUNFODJBMkI3RDE1QiIgc3RSZWY6ZG9jdW1lbnRJRD0ieG1w
LmRpZDoyMDBENzIxMzEyNUIxMUU4QTc5QUNFODJBMkI3RDE1QiIvPiA8L3JkZjpEZXNjcmlwdGlv
bj4gPC9yZGY6UkRGPiA8L3g6eG1wbWV0YT4gPD94cGFja2V0IGVuZD0iciI/PgH//v38+/r5+Pf2
9fTz8vHw7+7t7Ovq6ejn5uXk4+Lh4N/e3dzb2tnY19bV1NPS0dDPzs3My8rJyMfGxcTDwsHAv769
vLu6ubi3trW0s7KxsK+urayrqqmop6alpKOioaCfnp2cm5qZmJeWlZSTkpGQj46NjIuKiYiHhoWE
g4KBgH9+fXx7enl4d3Z1dHNycXBvbm1sa2ppaGdmZWRjYmFgX15dXFtaWVhXVlVUU1JRUE9OTUxL
SklIR0ZFRENCQUA/Pj08Ozo5ODc2NTQzMjEwLy4tLCsqKSgnJiUkIyIhIB8eHRwbGhkYFxYVFBMS
ERAPDg0MCwoJCAcGBQQDAgEAACH5BAEAABcALAAAAAAQABAAAAVn4CUiYlmS5qKaorqUURyb8kzd
OCXm+FX9wGCwBCkWBwNjkdVoNg4HZ4MlUlgDAauCKrJ4CQSvhcsoMwAAM4P1aLcLBXe7JKnbBQK7
/TLp+ycif34XDoWFJoaHIgmMVIwJLAZcF5IiIQA7');

my $cancel_image = $mw->Photo(-data => 'R0lGODlhEAAQAMQXAL8AANoAAKwAAKMAALcAAPHBwZkAAODCwtMAAKgAANsAAK8AAL4AALYAAKEA
AMUAAOkAAPgAAOIAAPEAAMwAAP8AAP///+7u7gAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACH/C1hN
UCBEYXRhWE1QPD94cGFja2V0IGJlZ2luPSLvu78iIGlkPSJXNU0wTXBDZWhpSHpyZVN6TlRjemtj
OWQiPz4gPHg6eG1wbWV0YSB4bWxuczp4PSJhZG9iZTpuczptZXRhLyIgeDp4bXB0az0iQWRvYmUg
WE1QIENvcmUgNS4wLWMwNjAgNjEuMTM0Nzc3LCAyMDEwLzAyLzEyLTE3OjMyOjAwICAgICAgICAi
PiA8cmRmOlJERiB4bWxuczpyZGY9Imh0dHA6Ly93d3cudzMub3JnLzE5OTkvMDIvMjItcmRmLXN5
bnRheC1ucyMiPiA8cmRmOkRlc2NyaXB0aW9uIHJkZjphYm91dD0iIiB4bWxuczp4bXA9Imh0dHA6
Ly9ucy5hZG9iZS5jb20veGFwLzEuMC8iIHhtbG5zOnhtcE1NPSJodHRwOi8vbnMuYWRvYmUuY29t
L3hhcC8xLjAvbW0vIiB4bWxuczpzdFJlZj0iaHR0cDovL25zLmFkb2JlLmNvbS94YXAvMS4wL3NU
eXBlL1Jlc291cmNlUmVmIyIgeG1wOkNyZWF0b3JUb29sPSJBZG9iZSBQaG90b3Nob3AgQ1M1IFdp
bmRvd3MiIHhtcE1NOkluc3RhbmNlSUQ9InhtcC5paWQ6RUU5MDE4OEMxQkRDMTFFODk1NkFGQzRC
NTYzQjg5MjQiIHhtcE1NOkRvY3VtZW50SUQ9InhtcC5kaWQ6RUU5MDE4OEQxQkRDMTFFODk1NkFG
QzRCNTYzQjg5MjQiPiA8eG1wTU06RGVyaXZlZEZyb20gc3RSZWY6aW5zdGFuY2VJRD0ieG1wLmlp
ZDpFRTkwMTg4QTFCREMxMUU4OTU2QUZDNEI1NjNCODkyNCIgc3RSZWY6ZG9jdW1lbnRJRD0ieG1w
LmRpZDpFRTkwMTg4QjFCREMxMUU4OTU2QUZDNEI1NjNCODkyNCIvPiA8L3JkZjpEZXNjcmlwdGlv
bj4gPC9yZGY6UkRGPiA8L3g6eG1wbWV0YT4gPD94cGFja2V0IGVuZD0iciI/PgH//v38+/r5+Pf2
9fTz8vHw7+7t7Ovq6ejn5uXk4+Lh4N/e3dzb2tnY19bV1NPS0dDPzs3My8rJyMfGxcTDwsHAv769
vLu6ubi3trW0s7KxsK+urayrqqmop6alpKOioaCfnp2cm5qZmJeWlZSTkpGQj46NjIuKiYiHhoWE
g4KBgH9+fXx7enl4d3Z1dHNycXBvbm1sa2ppaGdmZWRjYmFgX15dXFtaWVhXVlVUU1JRUE9OTUxL
SklIR0ZFRENCQUA/Pj08Ozo5ODc2NTQzMjEwLy4tLCsqKSgnJiUkIyIhIB8eHRwbGhkYFxYVFBMS
ERAPDg0MCwoJCAcGBQQDAgEAACH5BAEAABcALAAAAAAQABAAAAVi4CWOZCkWVaquaSFGcCzL4mTf
OC5CfO/7IolwSCSKFMhkIJBMihBQi3QqhSJElCyVmqWIHuAwABAOixjoNIGQTosa8LhAEI+LFvj8
YJDPixKAgQYGgYEiBw6JiouJByaPJCEAOw==');

my $processing_gif = $mw->Animation(-format => 'gif', -data => 'R0lGODlhHwAQAJEAABfQAKCgoO7u7gAAACH/C05FVFNDQVBFMi4wAwEAAAAh/wtYTVAgRGF0YVhN
UDw/eHBhY2tldCBiZWdpbj0i77u/IiBpZD0iVzVNME1wQ2VoaUh6cmVTek5UY3prYzlkIj8+IDx4
OnhtcG1ldGEgeG1sbnM6eD0iYWRvYmU6bnM6bWV0YS8iIHg6eG1wdGs9IkFkb2JlIFhNUCBDb3Jl
IDUuMC1jMDYwIDYxLjEzNDc3NywgMjAxMC8wMi8xMi0xNzozMjowMCAgICAgICAgIj4gPHJkZjpS
REYgeG1sbnM6cmRmPSJodHRwOi8vd3d3LnczLm9yZy8xOTk5LzAyLzIyLXJkZi1zeW50YXgtbnMj
Ij4gPHJkZjpEZXNjcmlwdGlvbiByZGY6YWJvdXQ9IiIgeG1sbnM6eG1wPSJodHRwOi8vbnMuYWRv
YmUuY29tL3hhcC8xLjAvIiB4bWxuczp4bXBNTT0iaHR0cDovL25zLmFkb2JlLmNvbS94YXAvMS4w
L21tLyIgeG1sbnM6c3RSZWY9Imh0dHA6Ly9ucy5hZG9iZS5jb20veGFwLzEuMC9zVHlwZS9SZXNv
dXJjZVJlZiMiIHhtcDpDcmVhdG9yVG9vbD0iQWRvYmUgUGhvdG9zaG9wIENTNSBXaW5kb3dzIiB4
bXBNTTpJbnN0YW5jZUlEPSJ4bXAuaWlkOjA2MEQ5RUM2MDkyNzExRTg5QUJEOEZBNEQwNDNGNTFB
IiB4bXBNTTpEb2N1bWVudElEPSJ4bXAuZGlkOjA2MEQ5RUM3MDkyNzExRTg5QUJEOEZBNEQwNDNG
NTFBIj4gPHhtcE1NOkRlcml2ZWRGcm9tIHN0UmVmOmluc3RhbmNlSUQ9InhtcC5paWQ6MDYwRDlF
QzQwOTI3MTFFODlBQkQ4RkE0RDA0M0Y1MUEiIHN0UmVmOmRvY3VtZW50SUQ9InhtcC5kaWQ6MDYw
RDlFQzUwOTI3MTFFODlBQkQ4RkE0RDA0M0Y1MUEiLz4gPC9yZGY6RGVzY3JpcHRpb24+IDwvcmRm
OlJERj4gPC94OnhtcG1ldGE+IDw/eHBhY2tldCBlbmQ9InIiPz4B//79/Pv6+fj39vX08/Lx8O/u
7ezr6uno5+bl5OPi4eDf3t3c29rZ2NfW1dTT0tHQz87NzMvKycjHxsXEw8LBwL++vby7urm4t7a1
tLOysbCvrq2sq6qpqKempaSjoqGgn56dnJuamZiXlpWUk5KRkI+OjYyLiomIh4aFhIOCgYB/fn18
e3p5eHd2dXRzcnFwb25tbGtqaWhnZmVkY2JhYF9eXVxbWllYV1ZVVFNSUVBPTk1MS0pJSEdGRURD
QkFAPz49PDs6OTg3NjU0MzIxMC8uLSwrKikoJyYlJCMiISAfHh0cGxoZGBcWFRQTEhEQDw4NDAsK
CQgHBgUEAwIBAAAh+QQEMgAAACwAAAAAHwAQAAACR5SPqYvhLh6UIUr7TNXx1N2Bn8dxYmeSJxqq
blqe8tvG9AjX6Jzj+sb7sXzA2/C4ChKRJGbPuTMqj9DlFDO5ULCfC+MLRhQAACH5BAQyAAAALAIA
AgADAAwAAAIFhI+pywUAIfkEBDIAAAAsCAACAAMADAAAAgWEj6nLBQAh+QQEMgAAACwOAAIAAwAM
AAACBYSPqcsFACH5BAQyAAAALBQAAgADAAwAAAIFhI+pywUAIfkEBDIAAAAsGgACAAMADAAAAgWE
j6nLBQAh+QQEMgAAACwCAAIAAwAMAAACBZSPqcsFACH5BAQyAAAALAgAAgADAAwAAAIFlI+pywUA
IfkEBDIAAAAsDgACAAMADAAAAgWUj6nLBQAh+QQEMgAAACwUAAIAAwAMAAACBZSPqcsFADs=');



#------------------------------------------#
#    CAPS markers filtration - window      #
#------------------------------------------#
my $Caps_Filtration_Window = $mw->Toplevel(-title => 'CAPS markers filtration');
$Caps_Filtration_Window->withdraw;
$Caps_Filtration_Window->iconimage($logo);
$Caps_Filtration_Window->resizable(1,1);
$Caps_Filtration_Window->geometry("+0+0");


#---------------------#
# NoteBook activation #
#---------------------#
my $nb = $Caps_Filtration_Window->NoteBook(
	-relief => 'raised',
	-bd => 1,
	-dynamicgeometry => 1
)->pack(-expand => 1, -padx => 2, -pady => 4, -fill => 'both');


#-------------------------#
# NoteBook pages creation #
#-------------------------#
my $genotype_filtration = $nb->add(
	'gf',
	-label => 'Filtration by genotype',
	-anchor => 'nw',
	-raisecmd => sub {
		$Caps_Filtration_Window->minsize(700,591);
		$Caps_Filtration_Window->maxsize(0,0);
		$mw->update;
	}
);

my $singleCut_filtration = $nb->add(
	'scf',
	-label => 'Filtration by cut-site',
	-anchor => 'nw',
	-raisecmd => sub {
		$Caps_Filtration_Window->minsize(420,150);
		$Caps_Filtration_Window->maxsize(420,150);
#		$Caps_Filtration_Window->maxsize(0,0);
		$mw->update;
	}
);

my $caps_to_fasta = $nb->add(
	'c2f',
	-label => 'CAPS to FASTA',
	-anchor => 'nw',
	-raisecmd => sub {
		$Caps_Filtration_Window->minsize(420,150);
		$Caps_Filtration_Window->maxsize(420,150);
#		$Caps_Filtration_Window->maxsize(0,0);
		$mw->update;
	}
);


#-------------------------------------#
# NoteBook - genotype filtration page #
#-------------------------------------#
my $cfw_gf_inputFile_check;
my $cfw_gf_inputFile_frame = $genotype_filtration->Frame->pack(-side => 'top', -fill => 'x');
my $cfw_gf_title_frame = $genotype_filtration->Frame->pack(-side => 'top', -fill => 'x');
my $cfw_gf_options_frame = $genotype_filtration->Frame->pack(-side => 'top', -fill => 'x');
my $cfw_gf_input_frame = $genotype_filtration->Frame->pack(-expand => 1, -side => 'top', -fill => 'both');

$cfw_gf_inputFile_frame->Label(-text => 'Please, select the VCF2CAPS output file with identified CAPS markers:')->pack(-side => 'top', -pady => 5, -padx => 5, -anchor => 'w');
$cfw_gf_inputFile_frame->Label(-text => 'VCF2CAPS output file:')->pack(-side => 'left', -padx => 5, -pady => 5);
my $cfw_gf_inputFile_entry = $cfw_gf_inputFile_frame->Entry(-insertwidth => 1, -width => 20,-textvariable => \$cfw_gf_input_file)->pack(-side => 'left');
my $cfw_gf_inputFile_chooseFile_button = $cfw_gf_inputFile_frame->Button(
	-text => 'D',
	-image => $folder_image,
	-command => sub { fileDialog($mw,"vcf2caps_output") }
)->pack(-side => 'left');
my $cfw_gf_inputFile_analyze_button;
my $cfw_gf_start_button;
$cfw_gf_inputFile_analyze_button = $cfw_gf_inputFile_frame->Button(
	-text => 'D',
	-image => $analyze_image,
	-command => sub {
		if (defined $cfw_gf_input_file and -e $cfw_gf_input_file and $jobID == 0)
		{
			start_vcf2capsOutput_check('gf');
			$cfw_gf_inputFile_analyze_button->configure(-state => 'disabled');
			$cfw_gf_inputFile_chooseFile_button->configure(-state => 'disabled');
		}
		elsif (defined $cfw_gf_input_file and !-e $cfw_gf_input_file)
		{
			$cfw_gf_inputFile_check->configure(-text => 'The file does not exist', -foreground => 'red', -image => '');
			$cfw_gf_start_button->configure(-state => 'disabled');
		}		
	}
)->pack(-side => 'left');
$cfw_gf_inputFile_check = $cfw_gf_inputFile_frame->Label(
	-text => 'none',
	-foreground => 'grey',
)->pack(-side => 'left', -padx => 5);


$cfw_gf_title_frame->Label(
	-text => 'To filter CAPS markers that differentiates specific individuals, paste into the following text fields space-, tab- or comma-separated names of individuals.

Below you can choose the way how CAPS markers will be selected:',
	-wraplength => 600,
	-justify => 'left'
	)->pack(-side => 'left', -padx => 5, -pady => 5);

my $cfw_gf_option_1 = $cfw_gf_options_frame->Radiobutton(
	-text => 'Only one group corresponds to the specific marker\'s genotype',
	-value => 1,
	-variable => \$oneFiltGroup_twoGenotypes
	)->pack(-side => 'left', -padx => 5, -pady => 5);
my $cfw_gf_option_2 = $cfw_gf_options_frame->Radiobutton(
	-text => 'Each group corresponds to the specific marker\'s genotype',
	-value => 0,
	-variable => \$oneFiltGroup_twoGenotypes
	)->pack(-side => 'left', -padx => 5, -pady => 5);
$cfw_gf_option_1->select();
	
my $cfw_gf_group_1_labelFrame = $cfw_gf_input_frame->Labelframe(-text => 'Group 1')->pack(-expand => 1, -fill => 'both', -padx => 5, -pady => 5);
my $cfw_gf_group_1_text = $cfw_gf_group_1_labelFrame->Scrolled('Text', -scrollbars => 'e', -insertwidth => 1, -wrap => 'word', -foreground => 'black' , -background => 'white', -relief => 'groove', -pady => 5, -padx => 5, -height => 5, width => 60)->pack(-expand => 1, -fill => 'both', -padx => 5, -pady => 5, -side => 'left');
$cfw_gf_group_1_text->bind('<KeyPress>', \&cfw_gf_groupsNo_checker);
my $cfw_gf_group1_properties_frame = $cfw_gf_group_1_labelFrame->Frame()->pack();
my $cfw_gf_group_1_maxError_entry = $cfw_gf_group1_properties_frame->Entry(
	-width => 5,
	-insertwidth => 2,
	-justify => 'right',
	-state => 'normal',
	-textvariable => \$cfw_gf_group_1_maxError_value
	)->pack(-padx => 2, -pady => 5, -side => 'left', -anchor => 'n');
$cfw_gf_group1_properties_frame->Label(-text => '%  Max mismatches')->pack(-pady => 5);;
$cfw_gf_group_1_labelFrame->Button(-text => 'Clear', -width => 4, -command => sub { $cfw_gf_group_1_text->delete('1.0', 'end') } )->pack(-padx => 2, -anchor => 'w');

my $cfw_gf_group_2_labelFrame = $cfw_gf_input_frame->Labelframe(-text => 'Group 2')->pack(-expand => 1, -fill => 'both', -padx => 5, -pady => 5);
my $cfw_gf_group_2_text = $cfw_gf_group_2_labelFrame->Scrolled('Text', -scrollbars => 'e', -insertwidth => 1, -wrap => 'word', -foreground => 'black' , -background => 'white', -relief => 'groove', -pady => 5, -padx => 5, -height => 5, width => 60)->pack(-expand => 1, -fill => 'both', -padx => 5, -pady => 5, -side => 'left');
$cfw_gf_group_2_text->bind('<KeyPress>', \&cfw_gf_groupsNo_checker);
my $cfw_gf_group2_properties_frame = $cfw_gf_group_2_labelFrame->Frame()->pack();
my $cfw_gf_group_2_maxError_entry = $cfw_gf_group2_properties_frame->Entry(
	-width => 5,
	-insertwidth => 2,
	-justify => 'right',
	-state => 'normal',
	-textvariable => \$cfw_gf_group_2_maxError_value
	)->pack(-padx => 2, -pady => 5, -side => 'left', -anchor => 'n');
$cfw_gf_group2_properties_frame->Label(-text => '%  Max mismatches')->pack(-pady => 5);
$cfw_gf_group_2_labelFrame->Button(-text => 'Clear', -width => 4, -command => sub { $cfw_gf_group_2_text->delete('1.0', 'end') } )->pack(-padx => 2, -anchor => 'w');


my $cfw_gf_group_3_labelFrame = $cfw_gf_input_frame->Labelframe(-text => 'Group 3')->pack(-expand => 1, -fill => 'both', -padx => 5, -pady => 5);
my $cfw_gf_group_3_text = $cfw_gf_group_3_labelFrame->Scrolled('Text', -scrollbars => 'e', -insertwidth => 1, -wrap => 'word', -foreground => 'black' , -background => 'white', -relief => 'groove', -pady => 5, -padx => 5, -height => 5, width => 60)->pack(-expand => 1, -fill => 'both', -padx => 5, -pady => 5, -side => 'left');

$cfw_gf_group_3_text->bind('<KeyPress>', \&cfw_gf_groupsNo_checker);
my $cfw_gf_group3_properties_frame = $cfw_gf_group_3_labelFrame->Frame()->pack();
my $cfw_gf_group_3_maxError_entry = $cfw_gf_group3_properties_frame->Entry(
	-width => 5,
	-insertwidth => 2,
	-justify => 'right',
	-state => 'normal',
	-textvariable => \$cfw_gf_group_3_maxError_value
	)->pack(-padx => 2, -pady => 5, -side => 'left', -anchor => 'n');
$cfw_gf_group3_properties_frame->Label(-text => '%  Max mismatches')->pack(-pady => 5);;
$cfw_gf_group_3_labelFrame->Button(-text => 'Clear', -width => 4, -command => sub { $cfw_gf_group_3_text->delete('1.0', 'end') } )->pack(-padx => 2, -anchor => 'w');

sub cfw_gf_groupsNo_checker
{
	my $cfw_group_1_value = $cfw_gf_group_1_text->get('1.0', 'end-1c');
	my $cfw_group_2_value = $cfw_gf_group_2_text->get('1.0', 'end-1c');
	my $cfw_group_3_value = $cfw_gf_group_3_text->get('1.0', 'end-1c');
	
	if ( $cfw_group_1_value ne "" and $cfw_group_2_value ne "" and $cfw_group_3_value ne "" )
	{
		$cfw_gf_option_2->configure(-value => 1);
		$cfw_gf_option_1->configure(-value => 0);
		$cfw_gf_option_1->configure(-state => 'disabled');
	}
	else
	{
		$cfw_gf_option_1->configure(-state => 'normal');
	}
}

my $cfw_gf_error_label = $cfw_gf_input_frame->Label();
my $cfw_gf_progress_frame = $cfw_gf_input_frame->Frame;
my $cfw_gf_stop_button = $cfw_gf_progress_frame->Button(-image => $cancel_image, -command => sub { $stop = 1 } )->pack(-side => 'left', -anchor => 'w');
my $cfw_gf_progress_textFrame = $cfw_gf_progress_frame->Text(-width => 13, -height => 1, -state => 'disabled')->pack(-side => 'left', -anchor => 'w', -padx => 5);
my $cfw_gf_result_label = $cfw_gf_progress_frame->Label()->pack(-side => 'left', -anchor => 'w');
my $cfw_gf_progressBar = $cfw_gf_progress_textFrame->ProgressBar(-variable => \$cfw_gf_CAPS_filtering_percent, -width => 14, -length => 90, -gap => 0, -from => 0, -to => 100, -foreground => 'blue', -troughcolor => 'white');
$cfw_gf_progress_textFrame->windowCreate('end', -window => $cfw_gf_progressBar);
my $cfw_gf_progress_label = $cfw_gf_progress_frame->Label()->pack(-side => 'left', -anchor => 'w');

$Caps_Filtration_Window->protocol('WM_DELETE_WINDOW' => sub { $Caps_Filtration_Window->withdraw } );

$cfw_gf_start_button = $cfw_gf_input_frame->Button(
	-text => 'Start filtration',
	-state => 'disabled',
	-command => sub {
		if ( $cfw_gf_group_1_maxError_value !~ /^[0-9]+\.?[0-9]*$/ or $cfw_gf_group_2_maxError_value !~ /^[0-9]+\.?[0-9]*$/ or $cfw_gf_group_3_maxError_value !~ /^[0-9]+\.?[0-9]*$/ )
		{
			curr_time();
			$terminal->insert('end', "Warning", 'warning');
			$terminal->insert('end', " - something is wrong with the marked parameter value. The parameter must be numerical.\n\n");
			$terminal->see('end');
			$cfw_gf_error_label->configure(-text => "Something is wrong with the marked parameter value. The parameter must be numerical.", -foreground => 'red');
			$cfw_gf_error_label->pack(-side => 'left', -padx => 5, -pady => 5, -anchor => 'w');
			
			my %group_number;
			$group_number{1}{value} = $cfw_gf_group_1_maxError_value;
			$group_number{1}{entry} = $cfw_gf_group_1_maxError_entry;
			$group_number{2}{value} = $cfw_gf_group_2_maxError_value;
			$group_number{2}{entry} = $cfw_gf_group_2_maxError_entry;
			$group_number{3}{value} = $cfw_gf_group_3_maxError_value;
			$group_number{3}{entry} = $cfw_gf_group_3_maxError_entry;
			
			foreach my $key ( keys (%group_number) )
			{
				if ( $group_number{$key}{value} !~ /^[0-9]+\.?[0-9]*$/ )
				{
					$group_number{$key}{entry}->configure(-background => 'red');
				}
				else
				{
					$group_number{$key}{entry}->configure(-background => 'white');
				}
			}			
			return;
		}
		else
		{
			$cfw_gf_group_1_maxError_entry->configure(-background => 'white');
			$cfw_gf_group_2_maxError_entry->configure(-background => 'white');
			$cfw_gf_group_3_maxError_entry->configure(-background => 'white');
			start_caps_filtration();
		}
	}
)->pack(-side => 'left', -padx => 5, -pady => 5, -anchor => 'w');



#---------------------------------------#
# NoteBook - single-cut filtration page #
#---------------------------------------#
my $cfw_scf_inputFile_check;
my $cfw_scf_inputFile_frame = $singleCut_filtration->Frame->pack(-side => 'top', -fill => 'x');
$cfw_scf_inputFile_frame->Label(-text => 'Please, select the VCF2CAPS output file with identified CAPS markers:')->pack(-side => 'top', -pady => 5, -padx => 5, -anchor => 'w');
$cfw_scf_inputFile_frame->Label(-text => 'VCF2CAPS output file:')->pack(-side => 'left', -padx => 5, -pady => 5);
my $cfw_scf_inputFile_entry = $cfw_scf_inputFile_frame->Entry(-insertwidth => 1, -width => 20,-textvariable => \$cfw_scf_input_file)->pack(-side => 'left');
my $cfw_scf_inputFile_chooseFile_button = $cfw_scf_inputFile_frame->Button(
	-text => 'D',
	-image => $folder_image,
	-command => sub { fileDialog($mw,"scf_vcf2caps_output") }
)->pack(-side => 'left');
my $cfw_scf_inputFile_analyze_button;
my $cfw_scf_start_button;

$cfw_scf_inputFile_analyze_button = $cfw_scf_inputFile_frame->Button(
	-text => 'D',
	-image => $analyze_image,
	-command => sub {
		if (defined $cfw_scf_input_file and -e $cfw_scf_input_file and $jobID == 0)
		{
			start_vcf2capsOutput_check('scf');
			$cfw_scf_inputFile_analyze_button->configure(-state => 'disabled');
			$cfw_scf_inputFile_chooseFile_button->configure(-state => 'disabled');
		}
		elsif (defined $cfw_scf_input_file and !-e $cfw_scf_input_file)
		{
			$cfw_scf_inputFile_check->configure(-text => 'The file does not exist', -foreground => 'red', -image => "");
			$cfw_scf_start_button->configure(-state => 'disabled');
		}
	}
)->pack(-side => 'left');

$cfw_scf_inputFile_check = $cfw_scf_inputFile_frame->Label(
	-wraplength => 120,
	-justify => 'left',
	-text => 'none',
	-foreground => 'grey',
)->pack(-side => 'left', -padx => 5);

my $cfw_scf_start_frame = $singleCut_filtration->Frame->pack(-side => 'top', -fill => 'x');

$cfw_scf_start_button = $cfw_scf_start_frame->Button(
	-text => 'Start filtration',
	-state => 'disabled',
	-command => sub {
		$cfw_scf_start_button->configure(-state => 'disabled');
		start_singleCut_filter();
	}
)->pack(-side => 'left', -padx => 5, -pady => 5, -anchor => 'w');

my $cfw_scf_error_label = $cfw_scf_start_frame->Label(-wraplength => 300, -justify => 'left');
my $cfw_scf_progress_frame = $cfw_scf_start_frame->Frame;
my $cfw_scf_stop_button = $cfw_scf_progress_frame->Button(-image => $cancel_image, -command => sub { $stop = 1 } )->pack(-side => 'left', -anchor => 'w');
my $cfw_scf_progress_textFrame = $cfw_scf_progress_frame->Text(-width => 13, -height => 1, -state => 'disabled')->pack(-side => 'left', -anchor => 'w', -padx => 5);
my $cfw_scf_result_label = $cfw_scf_progress_frame->Label()->pack(-side => 'left', -anchor => 'w');
my $cfw_scf_progressBar = $cfw_scf_progress_textFrame->ProgressBar(-variable => \$cfw_scf_CAPS_filtering_percent, -width => 14, -length => 90, -gap => 0, -from => 0, -to => 100, -foreground => 'blue', -troughcolor => 'white');
$cfw_scf_progress_textFrame->windowCreate('end', -window => $cfw_scf_progressBar);
my $cfw_scf_progress_label = $cfw_scf_progress_frame->Label()->pack(-side => 'left', -anchor => 'w');


#-------------------------------#
# NoteBook - CAPS to FASTA page #
#-------------------------------#
my $cfw_c2f_inputFile_check;
my $cfw_c2f_inputFile_frame = $caps_to_fasta->Frame->pack(-side => 'top', -fill => 'x');
$cfw_c2f_inputFile_frame->Label(-text => 'Please, select the VCF2CAPS output file with identified CAPS markers:')->pack(-side => 'top', -pady => 5, -padx => 5, -anchor => 'w');
$cfw_c2f_inputFile_frame->Label(-text => 'VCF2CAPS output file:')->pack(-side => 'left', -padx => 5, -pady => 5);
my $cfw_c2f_inputFile_entry = $cfw_c2f_inputFile_frame->Entry(-insertwidth => 1, -width => 20,-textvariable => \$cfw_c2f_input_file)->pack(-side => 'left');
my $cfw_c2f_inputFile_chooseFile_button = $cfw_c2f_inputFile_frame->Button(
	-text => 'D',
	-image => $folder_image,
	-command => sub { fileDialog($mw,"c2f_vcf2caps_output") }
)->pack(-side => 'left');
my $cfw_c2f_inputFile_analyze_button;
my $cfw_c2f_start_button;
$cfw_c2f_inputFile_analyze_button = $cfw_c2f_inputFile_frame->Button(
	-text => 'D',
	-image => $analyze_image,
	-command => sub {
		if (defined $cfw_c2f_input_file and -e $cfw_c2f_input_file and $jobID == 0)
		{
			start_vcf2capsOutput_check('c2f');
			$cfw_c2f_inputFile_analyze_button->configure(-state => 'disabled');
			$cfw_c2f_inputFile_chooseFile_button->configure(-state => 'disabled');
		}
		elsif (defined $cfw_c2f_input_file and !-e $cfw_c2f_input_file)
		{
			$cfw_c2f_inputFile_check->configure(-text => 'The file does not exist', -foreground => 'red', -image => "");
			$cfw_c2f_start_button->configure(-state => 'disabled');
		}
	}
)->pack(-side => 'left');
$cfw_c2f_inputFile_check = $cfw_c2f_inputFile_frame->Label(
	-wraplength => 120,
	-justify => 'left',
	-text => 'none',
	-foreground => 'grey',
)->pack(-side => 'left', -padx => 5);

my $cfw_c2f_start_frame = $caps_to_fasta->Frame->pack(-side => 'top', -fill => 'x');

$cfw_c2f_start_button = $cfw_c2f_start_frame->Button(
	-text => 'Start convertion',
	-state => 'disabled',
	-command => sub {
		$cfw_c2f_start_button->configure(-state => 'disabled');
		start_caps_to_fasta_convertion();
	}
	)->pack(-side => 'left', -padx => 5, -pady => 5, -anchor => 'w');

my $cfw_c2f_error_label = $cfw_c2f_start_frame->Label(-wraplength => 250, -justify => 'left');
my $cfw_c2f_progress_frame = $cfw_c2f_start_frame->Frame;
my $cfw_c2f_stop_button = $cfw_c2f_progress_frame->Button(-image => $cancel_image, -command => sub { $stop = 1 } )->pack(-side => 'left', -anchor => 'w');
my $cfw_c2f_progress_textFrame = $cfw_c2f_progress_frame->Text(-width => 13, -height => 1, -state => 'disabled')->pack(-side => 'left', -anchor => 'w', -padx => 5);
my $cfw_c2f_result_label = $cfw_c2f_progress_frame->Label()->pack(-side => 'left', -anchor => 'w');
my $cfw_c2f_progressBar = $cfw_c2f_progress_textFrame->ProgressBar(-variable => \$cfw_c2f_convertion_percent, -width => 14, -length => 90, -gap => 0, -from => 0, -to => 100, -foreground => 'blue', -troughcolor => 'white');
$cfw_c2f_progress_textFrame->windowCreate('end', -window => $cfw_c2f_progressBar);
my $cfw_c2f_progress_label = $cfw_c2f_progress_frame->Label()->pack(-side => 'left', -anchor => 'w');


#-----------------------------#
# 'About' window constructors #
#-----------------------------#
$mw->fontCreate('title', -family => 'arial', -size => 12, -weight => 'bold');
$mw->fontCreate('text', -family => 'arial', -size => 8);
$mw->fontCreate('text_b', -family => 'arial', -size => 8, -weight => 'bold');
$mw->fontCreate('hyper', -family => 'arial', -size => 8);

my $about_window = $mw->Toplevel(-title => 'About VCF2CAPS');
$about_window->withdraw;
$about_window->iconimage($logo);
$about_window->resizable(0,0);
my $about_text = $about_window->Text(-cursor => 'left_ptr', -width => 50, -height => 14, -insertwidth => 0, -insertontime => 0, -insertofftime => 0, -wrap => 'word', -background => 'gray95')->pack();
$about_text->tagConfigure('title_center', -justify => 'center', -font => 'title');
$about_text->tagConfigure('text_center', -justify => 'center', -font => 'text');
$about_text->tagConfigure('hyperlink', -underline => 0, -font => 'text', -foreground => 'blue', -justify => 'center');
$about_text->tagBind('hyperlink', '<Any-Enter>' => sub { $about_text->tagConfigure('hyperlink', -underline => 1, -font => 'text') } );
$about_text->tagBind('hyperlink', '<Any-Leave>' => sub { $about_text->tagConfigure('hyperlink', -underline => 0, -font => 'text') } );
$about_text->tagBind('hyperlink', '<Button-1>' => sub {
	open_hyperlink("https://github.com/Aviatore/VCF2CAPS.git");
} );
$about_text->insert('end',"\n");
$about_text->insert('end',"VCF2CAPS v2.0\n", 'title_center');
$about_text->insert('end',"Copyright \x{00A9} 2018 Wojciech Wesołowski\n\n", 'text_center');
$about_text->insert('end',"Free, open-source CAPS mining software from VCF files.\n", 'text_center');
$about_text->insert('end'," ", 'text_center');
$about_text->insert('end', "https://github.com/Aviatore/VCF2CAPS.git", 'hyperlink');
$about_text->insert('end',"\n\nUniversity of Agriculture in Krakow, Poland\n\n", 'text_center');
$about_text->insert('end',"Programmed by Wojciech Wesołowski\n", 'text_center');
$about_text->insert('end', "w.wesolowski\@protonmail.com\n\n", 'text_center');
$about_text->insert('end'," ", 'text_center');

my $button_OK = $about_text->Button(-text => 'OK', -width => 10, -command => sub { $about_window->withdraw } );
$about_text->windowCreate('text_center.last', -window => $button_OK, -align => 'center');
$about_window->protocol('WM_DELETE_WINDOW' => sub { $about_window->withdraw } );


#-------------------------------#
# 'License' window constructors #
#-------------------------------#
my $license_window = $mw->Toplevel(-title => 'License');
$license_window->withdraw;
$license_window->iconimage($logo);
my $license_text = $license_window->Scrolled('Text', -scrollbars => 'e', -padx => 10, -pady => 10, -cursor => 'left_ptr', -width => 100, -height => 28, -insertwidth => 0, -insertontime => 0, -insertofftime => 0, -wrap => 'word', -background => 'gray95')->pack(-fill => 'both', -expand => 1);
$license_text->tagConfigure('title_center', -justify => 'center', -font => 'title');
$license_text->tagConfigure('title_left', -font => 'title');
$license_text->tagConfigure('title_left_s', -font => 'text_b');
$license_text->tagConfigure('text', -font => 'text');
$license_text->tagConfigure('text_center', -font => 'text', -justify => 'center');
$license_text->insert('end',"\n");
$license_text->insert('end',"GNU GENERAL PUBLIC LICENSE\n", 'title_center');
$license_text->insert('end',"Version 3, 29 June 2007\n\n", 'text_center');
$license_text->insert('end',"Copyright \x{00A9} 2007 Free Software Foundation, Inc. <https://fsf.org/>\n", 'text');
$license_text->insert('end',"Everyone is permitted to copy and distribute verbatim copies of this license document, but changing it is not allowed.\n\n", 'text');
$license_text->insert('end',"Preamble\n\n", 'title_left');
$license_text->insert('end',"The GNU General Public License is a free, copyleft license for software and other kinds of works.

The licenses for most software and other practical works are designed to take away your freedom to share and change the works. By contrast, the GNU General Public License is intended to guarantee your freedom to share and change all versions of a program--to make sure it remains free software for all its users. We, the Free Software Foundation, use the GNU General Public License for most of our software; it applies also to any other work released this way by its authors. You can apply it to your programs, too.

When we speak of free software, we are referring to freedom, not price. Our General Public Licenses are designed to make sure that you have the freedom to distribute copies of free software (and charge for them if you wish), that you receive source code or can get it if you want it, that you can change the software or use pieces of it in new free programs, and that you know you can do these things.

To protect your rights, we need to prevent others from denying you these rights or asking you to surrender the rights. Therefore, you have certain responsibilities if you distribute copies of the software, or if you modify it: responsibilities to respect the freedom of others.

For example, if you distribute copies of such a program, whether gratis or for a fee, you must pass on to the recipients the same freedoms that you received. You must make sure that they, too, receive or can get the source code. And you must show them these terms so they know their rights.

Developers that use the GNU GPL protect your rights with two steps: (1) assert copyright on the software, and (2) offer you this License giving you legal permission to copy, distribute and/or modify it.

For the developers' and authors' protection, the GPL clearly explains that there is no warranty for this free software. For both users' and authors' sake, the GPL requires that modified versions be marked as changed, so that their problems will not be attributed erroneously to authors of previous versions.

Some devices are designed to deny users access to install or run modified versions of the software inside them, although the manufacturer can do so. This is fundamentally incompatible with the aim of protecting users' freedom to change the software. The systematic pattern of such abuse occurs in the area of products for individuals to use, which is precisely where it is most unacceptable. Therefore, we have designed this version of the GPL to prohibit the practice for those products. If such problems arise substantially in other domains, we stand ready to extend this provision to those domains in future versions of the GPL, as needed to protect the freedom of users.

Finally, every program is threatened constantly by software patents. States should not allow patents to restrict development and use of software on general-purpose computers, but in those that do, we wish to avoid the special danger that patents applied to a free program could make it effectively proprietary. To prevent this, the GPL assures that patents cannot be used to render the program non-free.

The precise terms and conditions for copying, distribution and modification follow.\n\n", 'text');
$license_text->insert('end',"TERMS AND CONDITIONS\n\n", 'title_left');
$license_text->insert('end',"0. Definitions.\n\n", 'title_left_s');
$license_text->insert('end',"\x{201F}This License\x{201D} refers to version 3 of the GNU General Public License.

\x{201F}Copyright\x{201D} also means copyright-like laws that apply to other kinds of works, such as semiconductor masks.

\x{201F}The Program\x{201D} refers to any copyrightable work licensed under this License. Each licensee is addressed as \x{201F}you\x{201D}. \x{201F}Licensees\x{201D} and \x{201F}recipients\x{201D} may be individuals or organizations.

To \x{201F}modify\x{201D} a work means to copy from or adapt all or part of the work in a fashion requiring copyright permission, other than the making of an exact copy. The resulting work is called a \x{201F}modified version\x{201D} of the earlier work or a work \x{201F}based on\x{201D} the earlier work.

A \x{201F}covered work\x{201D} means either the unmodified Program or a work based on the Program.

To \x{201F}propagate\x{201D} a work means to do anything with it that, without permission, would make you directly or secondarily liable for infringement under applicable copyright law, except executing it on a computer or modifying a private copy. Propagation includes copying, distribution (with or without modification), making available to the public, and in some countries other activities as well.

To \x{201F}convey\x{201D} a work means any kind of propagation that enables other parties to make or receive copies. Mere interaction with a user through a computer network, with no transfer of a copy, is not conveying.

An interactive user interface displays \x{201F}Appropriate Legal Notices\x{201D} to the extent that it includes a convenient and prominently visible feature that (1) displays an appropriate copyright notice, and (2) tells the user that there is no warranty for the work (except to the extent that warranties are provided), that licensees may convey the work under this License, and how to view a copy of this License. If the interface presents a list of user commands or options, such as a menu, a prominent item in the list meets this criterion.\n\n", 'text');
$license_text->insert('end',"1. Source Code.\n\n", 'title_left_s');
$license_text->insert('end',"The \x{201F}source code\x{201D} for a work means the preferred form of the work for making modifications to it. \x{201F}Object code\x{201D} means any non-source form of a work.

A \x{201F}Standard Interface\x{201D} means an interface that either is an official standard defined by a recognized standards body, or, in the case of interfaces specified for a particular programming language, one that is widely used among developers working in that language.

The \x{201F}System Libraries\x{201D} of an executable work include anything, other than the work as a whole, that (a) is included in the normal form of packaging a Major Component, but which is not part of that Major Component, and (b) serves only to enable use of the work with that Major Component, or to implement a Standard Interface for which an implementation is available to the public in source code form. A \x{201F}Major Component\x{201D}, in this context, means a major essential component (kernel, window system, and so on) of the specific operating system (if any) on which the executable work runs, or a compiler used to produce the work, or an object code interpreter used to run it.

The \x{201F}Corresponding Source\x{201D} for a work in object code form means all the source code needed to generate, install, and (for an executable work) run the object code and to modify the work, including scripts to control those activities. However, it does not include the work's System Libraries, or general-purpose tools or generally available free programs which are used unmodified in performing those activities but which are not part of the work. For example, Corresponding Source includes interface definition files associated with source files for the work, and the source code for shared libraries and dynamically linked subprograms that the work is specifically designed to require, such as by intimate data communication or control flow between those subprograms and other parts of the work.

The Corresponding Source need not include anything that users can regenerate automatically from other parts of the Corresponding Source.

The Corresponding Source for a work in source code form is that same work.
\n", 'text');
$license_text->insert('end',"2. Basic Permissions.\n\n", 'title_left_s');
$license_text->insert('end',"All rights granted under this License are granted for the term of copyright on the Program, and are irrevocable provided the stated conditions are met. This License explicitly affirms your unlimited permission to run the unmodified Program. The output from running a covered work is covered by this License only if the output, given its content, constitutes a covered work. This License acknowledges your rights of fair use or other equivalent, as provided by copyright law.

You may make, run and propagate covered works that you do not convey, without conditions so long as your license otherwise remains in force. You may convey covered works to others for the sole purpose of having them make modifications exclusively for you, or provide you with facilities for running those works, provided that you comply with the terms of this License in conveying all material for which you do not control copyright. Those thus making or running the covered works for you must do so exclusively on your behalf, under your direction and control, on terms that prohibit them from making any copies of your copyrighted material outside their relationship with you.

Conveying under any other circumstances is permitted solely under the conditions stated below. Sublicensing is not allowed; section 10 makes it unnecessary.
\n", 'text');
$license_text->insert('end',"3. Protecting Users' Legal Rights From Anti-Circumvention Law.\n\n", 'title_left_s');
$license_text->insert('end',"No covered work shall be deemed part of an effective technological measure under any applicable law fulfilling obligations under article 11 of the WIPO copyright treaty adopted on 20 December 1996, or similar laws prohibiting or restricting circumvention of such measures.

When you convey a covered work, you waive any legal power to forbid circumvention of technological measures to the extent such circumvention is effected by exercising rights under this License with respect to the covered work, and you disclaim any intention to limit operation or modification of the work as a means of enforcing, against the work's users, your or third parties' legal rights to forbid circumvention of technological measures.\n\n", 'text');
$license_text->insert('end',"4. Conveying Verbatim Copies.\n\n", 'title_left_s');
$license_text->insert('end',"You may convey verbatim copies of the Program's source code as you receive it, in any medium, provided that you conspicuously and appropriately publish on each copy an appropriate copyright notice; keep intact all notices stating that this License and any non-permissive terms added in accord with section 7 apply to the code; keep intact all notices of the absence of any warranty; and give all recipients a copy of this License along with the Program.

You may charge any price or no price for each copy that you convey, and you may offer support or warranty protection for a fee.\n\n", 'text');
$license_text->insert('end',"5. Conveying Modified Source Versions.\n\n", 'title_left_s');
$license_text->insert('end',"You may convey a work based on the Program, or the modifications to produce it from the Program, in the form of source code under the terms of section 4, provided that you also meet all of these conditions:

  a) The work must carry prominent notices stating that you modified it, and giving a relevant date.
  
  b) The work must carry prominent notices stating that it is released under this License and any conditions added under section 7. This requirement modifies the requirement in section 4 to “keep intact all notices”.
  
  c) You must license the entire work, as a whole, under this License to anyone who comes into possession of a copy. This License will therefore apply, along with any applicable section 7 additional terms, to the whole of the work, and all its parts, regardless of how they are packaged. This License gives no permission to license the work in any other way, but it does not invalidate such permission if you have separately received it.
  
  d) If the work has interactive user interfaces, each must display Appropriate Legal Notices; however, if the Program has interactive interfaces that do not display Appropriate Legal Notices, your work need not make them do so.
  
A compilation of a covered work with other separate and independent works, which are not by their nature extensions of the covered work, and which are not combined with it such as to form a larger program, in or on a volume of a storage or distribution medium, is called an “aggregate” if the compilation and its resulting copyright are not used to limit the access or legal rights of the compilation's users beyond what the individual works permit. Inclusion of a covered work in an aggregate does not cause this License to apply to the other parts of the aggregate.\n\n", 'text');
$license_text->insert('end',"6. Conveying Non-Source Forms.\n\n", 'title_left_s');
$license_text->insert('end',"You may convey a covered work in object code form under the terms of sections 4 and 5, provided that you also convey the machine-readable Corresponding Source under the terms of this License, in one of these ways:

  a) Convey the object code in, or embodied in, a physical product (including a physical distribution medium), accompanied by the Corresponding Source fixed on a durable physical medium customarily used for software interchange.
  
  b) Convey the object code in, or embodied in, a physical product (including a physical distribution medium), accompanied by a written offer, valid for at least three years and valid for as long as you offer spare parts or customer support for that product model, to give anyone who possesses the object code either (1) a copy of the Corresponding Source for all the software in the product that is covered by this License, on a durable physical medium customarily used for software interchange, for a price no more than your reasonable cost of physically performing this conveying of source, or (2) access to copy the Corresponding Source from a network server at no charge.
  
  c) Convey individual copies of the object code with a copy of the written offer to provide the Corresponding Source. This alternative is allowed only occasionally and noncommercially, and only if you received the object code with such an offer, in accord with subsection 6b.
  
  d) Convey the object code by offering access from a designated place (gratis or for a charge), and offer equivalent access to the Corresponding Source in the same way through the same place at no further charge. You need not require recipients to copy the Corresponding Source along with the object code. If the place to copy the object code is a network server, the Corresponding Source may be on a different server (operated by you or a third party) that supports equivalent copying facilities, provided you maintain clear directions next to the object code saying where to find the Corresponding Source. Regardless of what server hosts the Corresponding Source, you remain obligated to ensure that it is available for as long as needed to satisfy these requirements.
  
  e) Convey the object code using peer-to-peer transmission, provided you inform other peers where the object code and Corresponding Source of the work are being offered to the general public at no charge under subsection 6d.
  
A separable portion of the object code, whose source code is excluded from the Corresponding Source as a System Library, need not be included in conveying the object code work.

A \x{201F}User Product\x{201D} is either (1) a \x{201F}consumer product\x{201D}, which means any tangible personal property which is normally used for personal, family, or household purposes, or (2) anything designed or sold for incorporation into a dwelling. In determining whether a product is a consumer product, doubtful cases shall be resolved in favor of coverage. For a particular product received by a particular user, \x{201F}normally used\x{201D} refers to a typical or common use of that class of product, regardless of the status of the particular user or of the way in which the particular user actually uses, or expects or is expected to use, the product. A product is a consumer product regardless of whether the product has substantial commercial, industrial or non-consumer uses, unless such uses represent the only significant mode of use of the product.

\x{201F}Installation Information\x{201D} for a User Product means any methods, procedures, authorization keys, or other information required to install and execute modified versions of a covered work in that User Product from a modified version of its Corresponding Source. The information must suffice to ensure that the continued functioning of the modified object code is in no case prevented or interfered with solely because modification has been made.

If you convey an object code work under this section in, or with, or specifically for use in, a User Product, and the conveying occurs as part of a transaction in which the right of possession and use of the User Product is transferred to the recipient in perpetuity or for a fixed term (regardless of how the transaction is characterized), the Corresponding Source conveyed under this section must be accompanied by the Installation Information. But this requirement does not apply if neither you nor any third party retains the ability to install modified object code on the User Product (for example, the work has been installed in ROM).

The requirement to provide Installation Information does not include a requirement to continue to provide support service, warranty, or updates for a work that has been modified or installed by the recipient, or for the User Product in which it has been modified or installed. Access to a network may be denied when the modification itself materially and adversely affects the operation of the network or violates the rules and protocols for communication across the network.

Corresponding Source conveyed, and Installation Information provided, in accord with this section must be in a format that is publicly documented (and with an implementation available to the public in source code form), and must require no special password or key for unpacking, reading or copying.\n\n", 'text');
$license_text->insert('end',"7. Additional Terms.\n\n", 'title_left_s');
$license_text->insert('end',"\x{201F}Additional permissions\x{201D} are terms that supplement the terms of this License by making exceptions from one or more of its conditions. Additional permissions that are applicable to the entire Program shall be treated as though they were included in this License, to the extent that they are valid under applicable law. If additional permissions apply only to part of the Program, that part may be used separately under those permissions, but the entire Program remains governed by this License without regard to the additional permissions.

When you convey a copy of a covered work, you may at your option remove any additional permissions from that copy, or from any part of it. (Additional permissions may be written to require their own removal in certain cases when you modify the work.) You may place additional permissions on material, added by you to a covered work, for which you have or can give appropriate copyright permission.

Notwithstanding any other provision of this License, for material you add to a covered work, you may (if authorized by the copyright holders of that material) supplement the terms of this License with terms:

a) Disclaiming warranty or limiting liability differently from the terms of sections 15 and 16 of this License; or

b) Requiring preservation of specified reasonable legal notices or author attributions in that material or in the Appropriate Legal Notices displayed by works containing it; or

c) Prohibiting misrepresentation of the origin of that material, or requiring that modified versions of such material be marked in reasonable ways as different from the original version; or

d) Limiting the use for publicity purposes of names of licensors or authors of the material; or

e) Declining to grant rights under trademark law for use of some trade names, trademarks, or service marks; or

f) Requiring indemnification of licensors and authors of that material by anyone who conveys the material (or modified versions of it) with contractual assumptions of liability to the recipient, for any liability that these contractual assumptions directly impose on those licensors and authors.

All other non-permissive additional terms are considered \x{201F}further restrictions\x{201D} within the meaning of section 10. If the Program as you received it, or any part of it, contains a notice stating that it is governed by this License along with a term that is a further restriction, you may remove that term. If a license document contains a further restriction but permits relicensing or conveying under this License, you may add to a covered work material governed by the terms of that license document, provided that the further restriction does not survive such relicensing or conveying.

If you add terms to a covered work in accord with this section, you must place, in the relevant source files, a statement of the additional terms that apply to those files, or a notice indicating where to find the applicable terms.

Additional terms, permissive or non-permissive, may be stated in the form of a separately written license, or stated as exceptions; the above requirements apply either way.\n\n", 'text');
$license_text->insert('end',"8. Termination.\n\n", 'title_left_s');
$license_text->insert('end',"You may not propagate or modify a covered work except as expressly provided under this License. Any attempt otherwise to propagate or modify it is void, and will automatically terminate your rights under this License (including any patent licenses granted under the third paragraph of section 11).

However, if you cease all violation of this License, then your license from a particular copyright holder is reinstated (a) provisionally, unless and until the copyright holder explicitly and finally terminates your license, and (b) permanently, if the copyright holder fails to notify you of the violation by some reasonable means prior to 60 days after the cessation.

Moreover, your license from a particular copyright holder is reinstated permanently if the copyright holder notifies you of the violation by some reasonable means, this is the first time you have received notice of violation of this License (for any work) from that copyright holder, and you cure the violation prior to 30 days after your receipt of the notice.

Termination of your rights under this section does not terminate the licenses of parties who have received copies or rights from you under this License. If your rights have been terminated and not permanently reinstated, you do not qualify to receive new licenses for the same material under section 10.\n\n", 'text');
$license_text->insert('end',"9. Acceptance Not Required for Having Copies.\n\n", 'title_left_s');
$license_text->insert('end',"You are not required to accept this License in order to receive or run a copy of the Program. Ancillary propagation of a covered work occurring solely as a consequence of using peer-to-peer transmission to receive a copy likewise does not require acceptance. However, nothing other than this License grants you permission to propagate or modify any covered work. These actions infringe copyright if you do not accept this License. Therefore, by modifying or propagating a covered work, you indicate your acceptance of this License to do so.\n\n", 'text');
$license_text->insert('end',"10. Automatic Licensing of Downstream Recipients.\n\n", 'title_left_s');
$license_text->insert('end',"Each time you convey a covered work, the recipient automatically receives a license from the original licensors, to run, modify and propagate that work, subject to this License. You are not responsible for enforcing compliance by third parties with this License.

An \x{201F}entity transaction\x{201D} is a transaction transferring control of an organization, or substantially all assets of one, or subdividing an organization, or merging organizations. If propagation of a covered work results from an entity transaction, each party to that transaction who receives a copy of the work also receives whatever licenses to the work the party's predecessor in interest had or could give under the previous paragraph, plus a right to possession of the Corresponding Source of the work from the predecessor in interest, if the predecessor has it or can get it with reasonable efforts.

You may not impose any further restrictions on the exercise of the rights granted or affirmed under this License. For example, you may not impose a license fee, royalty, or other charge for exercise of rights granted under this License, and you may not initiate litigation (including a cross-claim or counterclaim in a lawsuit) alleging that any patent claim is infringed by making, using, selling, offering for sale, or importing the Program or any portion of it.\n\n", 'text');
$license_text->insert('end',"11. Patents.\n\n", 'title_left_s');
$license_text->insert('end',"A \x{201F}contributor\x{201D} is a copyright holder who authorizes use under this License of the Program or a work on which the Program is based. The work thus licensed is called the contributor's \x{201F}contributor version\x{201D}.

A contributor's \x{201F}essential patent claims\x{201D} are all patent claims owned or controlled by the contributor, whether already acquired or hereafter acquired, that would be infringed by some manner, permitted by this License, of making, using, or selling its contributor version, but do not include claims that would be infringed only as a consequence of further modification of the contributor version. For purposes of this definition, \x{201F}control\x{201D} includes the right to grant patent sublicenses in a manner consistent with the requirements of this License.

Each contributor grants you a non-exclusive, worldwide, royalty-free patent license under the contributor's essential patent claims, to make, use, sell, offer for sale, import and otherwise run, modify and propagate the contents of its contributor version.

In the following three paragraphs, a \x{201F}patent license\x{201D} is any express agreement or commitment, however denominated, not to enforce a patent (such as an express permission to practice a patent or covenant not to sue for patent infringement). To \x{201F}grant\x{201D} such a patent license to a party means to make such an agreement or commitment not to enforce a patent against the party.

If you convey a covered work, knowingly relying on a patent license, and the Corresponding Source of the work is not available for anyone to copy, free of charge and under the terms of this License, through a publicly available network server or other readily accessible means, then you must either (1) cause the Corresponding Source to be so available, or (2) arrange to deprive yourself of the benefit of the patent license for this particular work, or (3) arrange, in a manner consistent with the requirements of this License, to extend the patent license to downstream recipients. \x{201F}Knowingly relying\x{201D} means you have actual knowledge that, but for the patent license, your conveying the covered work in a country, or your recipient's use of the covered work in a country, would infringe one or more identifiable patents in that country that you have reason to believe are valid.

If, pursuant to or in connection with a single transaction or arrangement, you convey, or propagate by procuring conveyance of, a covered work, and grant a patent license to some of the parties receiving the covered work authorizing them to use, propagate, modify or convey a specific copy of the covered work, then the patent license you grant is automatically extended to all recipients of the covered work and works based on it.

A patent license is \x{201F}discriminatory\x{201D} if it does not include within the scope of its coverage, prohibits the exercise of, or is conditioned on the non-exercise of one or more of the rights that are specifically granted under this License. You may not convey a covered work if you are a party to an arrangement with a third party that is in the business of distributing software, under which you make payment to the third party based on the extent of your activity of conveying the work, and under which the third party grants, to any of the parties who would receive the covered work from you, a discriminatory patent license (a) in connection with copies of the covered work conveyed by you (or copies made from those copies), or (b) primarily for and in connection with specific products or compilations that contain the covered work, unless you entered into that arrangement, or that patent license was granted, prior to 28 March 2007.

Nothing in this License shall be construed as excluding or limiting any implied license or other defenses to infringement that may otherwise be available to you under applicable patent law.\n\n", 'text');
$license_text->insert('end',"12. No Surrender of Others' Freedom.\n\n", 'title_left_s');
$license_text->insert('end',"If conditions are imposed on you (whether by court order, agreement or otherwise) that contradict the conditions of this License, they do not excuse you from the conditions of this License. If you cannot convey a covered work so as to satisfy simultaneously your obligations under this License and any other pertinent obligations, then as a consequence you may not convey it at all. For example, if you agree to terms that obligate you to collect a royalty for further conveying from those to whom you convey the Program, the only way you could satisfy both those terms and this License would be to refrain entirely from conveying the Program.\n\n", 'text');
$license_text->insert('end',"13. Use with the GNU Affero General Public License.\n\n", 'title_left_s');
$license_text->insert('end',"Notwithstanding any other provision of this License, you have permission to link or combine any covered work with a work licensed under version 3 of the GNU Affero General Public License into a single combined work, and to convey the resulting work. The terms of this License will continue to apply to the part which is the covered work, but the special requirements of the GNU Affero General Public License, section 13, concerning interaction through a network will apply to the combination as such.\n\n", 'text');
$license_text->insert('end',"14. Revised Versions of this License.\n\n", 'title_left_s');
$license_text->insert('end',"The Free Software Foundation may publish revised and/or new versions of the GNU General Public License from time to time. Such new versions will be similar in spirit to the present version, but may differ in detail to address new problems or concerns.

Each version is given a distinguishing version number. If the Program specifies that a certain numbered version of the GNU General Public License “or any later version” applies to it, you have the option of following the terms and conditions either of that numbered version or of any later version published by the Free Software Foundation. If the Program does not specify a version number of the GNU General Public License, you may choose any version ever published by the Free Software Foundation.

If the Program specifies that a proxy can decide which future versions of the GNU General Public License can be used, that proxy's public statement of acceptance of a version permanently authorizes you to choose that version for the Program.

Later license versions may give you additional or different permissions. However, no additional obligations are imposed on any author or copyright holder as a result of your choosing to follow a later version.\n\n", 'text');
$license_text->insert('end',"15. Disclaimer of Warranty.\n\n", 'title_left_s');
$license_text->insert('end',"THERE IS NO WARRANTY FOR THE PROGRAM, TO THE EXTENT PERMITTED BY APPLICABLE LAW. EXCEPT WHEN OTHERWISE STATED IN WRITING THE COPYRIGHT HOLDERS AND/OR OTHER PARTIES PROVIDE THE PROGRAM “AS IS” WITHOUT WARRANTY OF ANY KIND, EITHER EXPRESSED OR IMPLIED, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE. THE ENTIRE RISK AS TO THE QUALITY AND PERFORMANCE OF THE PROGRAM IS WITH YOU. SHOULD THE PROGRAM PROVE DEFECTIVE, YOU ASSUME THE COST OF ALL NECESSARY SERVICING, REPAIR OR CORRECTION.\n\n", 'text');
$license_text->insert('end',"16. Limitation of Liability.\n\n", 'title_left_s');
$license_text->insert('end',"IN NO EVENT UNLESS REQUIRED BY APPLICABLE LAW OR AGREED TO IN WRITING WILL ANY COPYRIGHT HOLDER, OR ANY OTHER PARTY WHO MODIFIES AND/OR CONVEYS THE PROGRAM AS PERMITTED ABOVE, BE LIABLE TO YOU FOR DAMAGES, INCLUDING ANY GENERAL, SPECIAL, INCIDENTAL OR CONSEQUENTIAL DAMAGES ARISING OUT OF THE USE OR INABILITY TO USE THE PROGRAM (INCLUDING BUT NOT LIMITED TO LOSS OF DATA OR DATA BEING RENDERED INACCURATE OR LOSSES SUSTAINED BY YOU OR THIRD PARTIES OR A FAILURE OF THE PROGRAM TO OPERATE WITH ANY OTHER PROGRAMS), EVEN IF SUCH HOLDER OR OTHER PARTY HAS BEEN ADVISED OF THE POSSIBILITY OF SUCH DAMAGES.\n\n", 'text');
$license_text->insert('end',"17. Interpretation of Sections 15 and 16.\n\n", 'title_left_s');
$license_text->insert('end',"If the disclaimer of warranty and limitation of liability provided above cannot be given local legal effect according to their terms, reviewing courts shall apply local law that most closely approximates an absolute waiver of all civil liability in connection with the Program, unless a warranty or assumption of liability accompanies a copy of the Program in return for a fee.

END OF TERMS AND CONDITIONS\n\n", 'text');

$license_text->insert('end'," ", 'text_center');

my $button_OK_lic = $license_text->Button(-text => 'Close', -width => 10, -command => sub { $license_window->withdraw } );
$license_text->windowCreate('text_center.last', -window => $button_OK_lic, -align => 'center');
$license_window->protocol('WM_DELETE_WINDOW' => sub { $license_window->withdraw } );


#------------#
# Menu START #
#------------#
$mw->configure(-menu => my $menubar = $mw->Menu);
my $file_menu = $menubar->cascade(-label => '~File');
my $help_menu = $menubar->cascade(-label => '~Help');

my $new_working_directory = $file_menu->command(-label => 'New working directory', -underline => 0, -command => sub { new_working_directory(\$working_dir) } );
my $download_db = $file_menu->command(-label => 'Download enzyme database', -underline => 0, -command => \&start_download_enzyme_db);
my $save_selected_enzymes = $file_menu->command(-label => 'Save selected enzymes', -underline => 0, -command => sub { fileDialog_save_enzyme($mw) } );
my $load_selected_enzymes = $file_menu->command(-label => 'Load selected enzymes', -underline => 0, -command => sub { fileDialog_load_enzyme($mw) });

$file_menu->separator;
$file_menu->command(-label => 'Exit', -underline => 0, -command => sub { exit });

my $license = $help_menu->command(-label => 'License information', -underline => 0, -command => sub { $license_text->yview(moveto => 0), $license_window->deiconify; $license_window->raise });
my $about = $help_menu->command(-label => 'About VCF2CAPS', -underline => 0, -command => sub { $about_window->deiconify; $about_window->raise } );


#-------------------------#
# Main window constructor #
#-------------------------#
my $top = $mw->Frame->pack(-side => 'top', -fill => 'x');
my $bottom = $mw->Frame->pack(-side => 'top', -expand => 1, -fill => 'both');

$bottom->Label(-text => 'Running information:')->pack(-side => 'top', -anchor => 'w', -padx => 5);

$terminal = $bottom->Scrolled('Text', -scrollbars => 'e', -insertwidth => 0, -insertontime => 0, -insertofftime => 0, -wrap => 'word', -foreground => 'gray95' , -background => 'black', -relief => 'groove', -pady => 5, -padx => 5, -height => 20);

$terminal->tagConfigure('warning', -foreground => 'red');
$terminal->tagConfigure('warning_p', -foreground => 'red', -spacing3 => 0);
$terminal->tagConfigure('mark', -foreground => 'yellow');
$terminal->tagConfigure('date', -foreground => 'SlateGrey');
$terminal->tagConfigure('paragraph', -spacing3 => 0);
$terminal->tagConfigure('stickUP', -spacing1 => 0);
$terminal->tagConfigure('percent');

my $L_frame = $top->Frame->pack(
	-side => 'left',
	-anchor => 'w'
);
my $R_frame = $top->Frame->pack(
	-side => 'left',
	-padx => 20
);

my $R_frame_allEnzymes_frame = $R_frame->Frame->pack(-side => 'left');
my $R_frame_buttons_frame = $R_frame->Frame->pack(-side => 'left');
my $R_frame_selEnzymes_frame = $R_frame->Frame->pack(-side => 'left');
$R_frame_allEnzymes_frame->Label(-text => 'All enzymes')->pack(-side => 'top',  -anchor => 'w', -padx => 10);
$R_frame_selEnzymes_frame->Label(-text => 'Selected enzymes')->pack(-side => 'top');
my $R_frame_selEnzymes_listBox;

my $R_frame_allEnzymes_listBox = $R_frame_allEnzymes_frame->Scrolled(
	'Listbox',
	-scrollbars => 'e',
	-height => 10,
	-width => 12,
	-selectmode => 'extended',
	-exportselection => 1,
)->pack(-side => 'top');
my $R_frame_allEnzymes_No = $R_frame_allEnzymes_frame->Label(-text => '0 enzymes')->pack(-side => 'top', -anchor => 'w', -padx => 10);
my $R_frame_selEnzymes_No;
my $R_frame_select_Button = $R_frame_buttons_frame->Button(
	-text => '>',
	-width => 2,
	-command => sub {
		my @selected_enz_index = $R_frame_allEnzymes_listBox->curselection();
		foreach my $index (@selected_enz_index)
		{
			push @selected_enz_names, $R_frame_allEnzymes_listBox->get($index);
		}
		@selected_enz_names = uniq(@selected_enz_names);
		@selected_enz_names = sort{$a cmp $b} @selected_enz_names;
		$R_frame_selEnzymes_listBox->delete(0, 'end');
		foreach my $enzyme_name (@selected_enz_names)
		{
			$R_frame_selEnzymes_listBox->insert('end', $enzyme_name);
		}
		
		$R_frame_selEnzymes_No->configure(-text => scalar(@selected_enz_names) . " enzymes" );
	}
)->pack(
	-anchor => 'c'
);
my $R_frame_deselect_Button = $R_frame_buttons_frame->Button(
	-text => '<',
	-width => 2,
	-command => sub {
		my @selected_enz_index = $R_frame_selEnzymes_listBox->curselection();
		foreach my $index (@selected_enz_index)
		{
			@selected_enz_names = grep { $_ ne $R_frame_selEnzymes_listBox->get($index) } @selected_enz_names;
		}
		$R_frame_selEnzymes_listBox->delete(0, 'end');
		foreach my $enzyme_name (@selected_enz_names)
		{
			$R_frame_selEnzymes_listBox->insert('end', $enzyme_name);
		}
		
		$R_frame_selEnzymes_No->configure(-text => scalar(@selected_enz_names) . " enzymes" );
	}
)->pack(
	-anchor => 'c'
);
my $R_frame_selectAll_Button = $R_frame_buttons_frame->Button(
	-text => '>>',
	-width => 2,
	-command => sub {
		my $all_enz_index_size = $R_frame_allEnzymes_listBox->size();
		for (my $i = 0; $i < $all_enz_index_size; $i++)
		{
			push @selected_enz_names, $R_frame_allEnzymes_listBox->get($i);
		}
		@selected_enz_names = uniq(@selected_enz_names);
		@selected_enz_names = sort{$a cmp $b} @selected_enz_names;
		$R_frame_selEnzymes_listBox->delete(0, 'end');
		foreach my $enzyme_name (@selected_enz_names)
		{
			$R_frame_selEnzymes_listBox->insert('end', $enzyme_name);
		}
		
		$R_frame_selEnzymes_No->configure(-text => scalar(@selected_enz_names) . " enzymes" );
	}
)->pack(
	-anchor => 'c'
);
my $R_frame_deselectAll_Button = $R_frame_buttons_frame->Button(
	-text => '<<',
	-width => 2,
	-command => sub {
		$R_frame_selEnzymes_listBox->delete(0, 'end');
		@selected_enz_names = ();
		
		$R_frame_selEnzymes_No->configure(-text => scalar(@selected_enz_names) . " enzymes" );
	}
)->pack(
	-anchor => 'c'
);
$R_frame_selEnzymes_listBox = $R_frame_selEnzymes_frame->Scrolled(
	'Listbox',
	-scrollbars => 'e',
	-height => 10,
	-width => 12,
	-exportselection => 1
)->pack(-side => 'top');

$R_frame_selEnzymes_No = $R_frame_selEnzymes_frame->Label(-text => '0 enzymes')->pack(-side => 'top', -anchor => 'w', -padx => 10);

my $R_frame_enzymeProperties_label_listBox;
my $R_frame_enzymeProperties_values_text;
my $R_frame_enzymeProperties_label_text;
my $enzyme_description = 0;
$R_frame_selEnzymes_listBox->bind('<<ListboxSelect>>' => sub {
	my $selected_enz_index = $R_frame_selEnzymes_listBox->curselection();

	if (defined $selected_enz_index)
	{
		$selected_enz_name = $R_frame_selEnzymes_listBox->get($selected_enz_index);
	}
	
	if ($enzyme_description == 0 and defined $selected_enz_index)
	{
		$R_frame_enzymeProperties_label_text->pack(-side => 'left', -fill => 'x');
		$R_frame_enzymeProperties_values_text->pack(-side => 'left', -fill => 'x');
		$enzyme_description = 1;
	}
	
	if ($comp_state == 1 and $iso_state == 0 and defined $selected_enz_index)
	{
		$R_frame_enzymeProperties_values_text->delete('0.0', 'end');
		
		if (( split("\t", $enzymes_db{$selected_enz_name}) )[4] eq "null")
		{
			$R_frame_enzymeProperties_values_text->insert('end', "None");
			$R_frame_enzymeProperties_label_text->tagConfigure('comp', -foreground => 'red');
			$R_frame_enzymeProperties_label_text->tagConfigure('iso', -foreground => 'black');
		}
		else
		{
			my $first = 0;
			foreach my $company_ID ( split("", ( split("\t", $enzymes_db{$selected_enz_name}) )[4] ) )
			{
				F: foreach my $company ( split("\t", $enzymes_db{companies}) )
				{
					if ( ( split(",", $company) )[0] eq  $company_ID)
					{
						my $company_name = ( split(",", $company) )[1];
						if ($first == 0) { $R_frame_enzymeProperties_values_text->insert('end', "$company_name") }
						else { $R_frame_enzymeProperties_values_text->insert('end', "\n$company_name") }
						$first++;
						last F;
					}
				}
				$R_frame_enzymeProperties_label_text->tagConfigure('comp', -foreground => 'red');
				$R_frame_enzymeProperties_label_text->tagConfigure('iso', -foreground => 'black');
			}
		}
	}
	elsif ($comp_state == 0 and $iso_state == 1 and defined $selected_enz_index)
	{
		$R_frame_enzymeProperties_values_text->delete('0.0', 'end');
		
		if (( split("\t", $enzymes_db{$selected_enz_name}) )[3] eq "null")
		{
			$R_frame_enzymeProperties_values_text->insert('end', "None");
			$R_frame_enzymeProperties_label_text->tagConfigure('iso', -foreground => 'red');
			$R_frame_enzymeProperties_label_text->tagConfigure('comp', -foreground => 'black');
		}
		else
		{
			my $first = 0;
			foreach my $isoschisomer ( split(",", ( split("\t", $enzymes_db{$selected_enz_name}) )[3] ) )
			{
				if ($first == 0) { $R_frame_enzymeProperties_values_text->insert('end', "$isoschisomer") }
				else { $R_frame_enzymeProperties_values_text->insert('end', "\n$isoschisomer") }
				$first++;
				
				$R_frame_enzymeProperties_label_text->tagConfigure('iso', -foreground => 'red');
				$R_frame_enzymeProperties_label_text->tagConfigure('comp', -foreground => 'black');
			}
		}
	}
	elsif ($comp_state == 0 and $iso_state == 0 and defined $selected_enz_index)
	{
		$R_frame_enzymeProperties_values_text->delete('0.0', 'end');
		
		my @data = split("\t", $enzymes_db{$selected_enz_name});
		my @seq;
		for ( my $i = 0; $i < scalar( split("", $data[1]) ); $i++ )
		{
			if ($data[0] == $i) { push @seq, ("'", ( split("", $data[1]) )[$i] ) }
			else
			{
				push @seq, ( split("", $data[1]) )[$i]
			}
		}
		$R_frame_enzymeProperties_values_text->delete('0.0', 'end');
		$R_frame_enzymeProperties_values_text->insert('end', $selected_enz_name . "\n");
		$R_frame_enzymeProperties_values_text->insert('end', $data[-1] . "\n");
		$R_frame_enzymeProperties_values_text->insert('end', $data[0] . "\n");
		$R_frame_enzymeProperties_values_text->insert('end', $data[2]);
		
		$R_frame_enzymeProperties_label_text->tagConfigure('comp', -foreground => 'black');
	}
});

my $R_frame_parametersContainer_frame = $R_frame->Frame->pack(-side => 'left', -fill => 'x');
my $R_frame_parameters_frame = $R_frame_parametersContainer_frame->Frame->pack(-side => 'top', -fill => 'both');
my $R_frame_parameters_commercially_checkButton = $R_frame_parameters_frame->Checkbutton(
	-text => 'Only commercially available enzymes',
	-variable => \$onlyCommercially
)->pack(-side => 'top', -anchor => 'w');
my $R_frame_parameters_4cutters_checkButton = $R_frame_parameters_frame->Checkbutton(
	-text => '4-base cutters',
	-variable => \$cutters4
)->pack(-side => 'top', -anchor => 'w');
my $R_frame_parameters_5cutters_checkButton = $R_frame_parameters_frame->Checkbutton(
	-text => '5-base cutters',
	-variable => \$cutters5
)->pack(-side => 'top', -anchor => 'w');
my $R_frame_parameters_6cutters_checkButton = $R_frame_parameters_frame->Checkbutton(
	-text => '6-base cutters',
	-variable => \$cutters6
)->pack(-side => 'top', -anchor => 'w');
my $R_frame_parameters_custom_entry;
my $R_frame_parameters_custom_frame = $R_frame_parameters_frame->Frame->pack(-side => 'top', -anchor => 'w');
my $R_frame_parameters_custom_checkButton = $R_frame_parameters_custom_frame->Checkbutton(
	-text => 'Custom',
	-variable => \$custom,
	-command => sub {
		if ($custom == 1)
		{
			$R_frame_parameters_4cutters_checkButton->configure(-state => 'disabled');
			$R_frame_parameters_5cutters_checkButton->configure(-state => 'disabled');
			$R_frame_parameters_6cutters_checkButton->configure(-state => 'disabled');
			$R_frame_parameters_custom_entry->configure(-state => 'normal');
			$custom_value = "";			
		}
		else
		{
			$R_frame_parameters_4cutters_checkButton->configure(-state => 'normal');
			$R_frame_parameters_5cutters_checkButton->configure(-state => 'normal');
			$R_frame_parameters_6cutters_checkButton->configure(-state => 'normal');
			$R_frame_parameters_custom_entry->configure(-state => 'disabled');
			$custom_value = "";
		}
	}
)->pack(-side => 'left', -anchor => 'w');
$R_frame_parameters_custom_entry = $R_frame_parameters_custom_frame->Entry(
	-width => 5,
	-insertwidth => 2,
	-justify => 'right',
	-textvariable => \$custom_value,
	-state => 'disabled'
)->pack(-side => 'left');
$R_frame_parameters_custom_frame->Label(
	-text => '-base cutters',
	-foreground => 'black'
)->pack(-side => 'left');
my $R_frame_parameters_applyFilters_Button = $R_frame_parameters_frame->Button(
	-text => 'Apply filters',
	-command => sub {
		if (scalar(@selected_enz_names) > 0 )
		{
			my @filtered_enzymes;
			my @cutters_type;
			
			if ($custom == 1 and $custom_value ne "")
			{
				foreach my $cutter ( split(",", $custom_value) )
				{
					if ($cutter =~ /^[0-9]+$/)
					{
						push @cutters_type, $cutter;
					}
				}
			}
			elsif ($cutters4 == 1 or $cutters5 == 1 or $cutters6 == 1)
			{
				if ($cutters4 == 1) { push @cutters_type, "4" }
				if ($cutters5 == 1) { push @cutters_type, "5" }
				if ($cutters6 == 1) { push @cutters_type, "6" }
			}
		
			if ( scalar(@cutters_type) > 0 )
			{
				foreach my $enzyme_name (@selected_enz_names)
				{
					my $sequence = ( split("\t", $enzymes_db{$enzyme_name}) )[1];
					foreach my $cutter (@cutters_type)
					{
						if ( scalar( split("", $sequence) ) == $cutter )
						{
							if ($onlyCommercially == 1)
							{
								if ( ( split("\t", $enzymes_db{$enzyme_name}) )[4] ne "null" )
								{
									push @filtered_enzymes, $enzyme_name;
								}
							}
							else
							{
								push @filtered_enzymes, $enzyme_name;
							}
						}
					}
				}
			}
			elsif ($onlyCommercially == 1)
			{
				foreach my $enzyme_name (@selected_enz_names)
				{
					if ( ( split("\t", $enzymes_db{$enzyme_name}) )[4] ne "null" )
					{
						push @filtered_enzymes, $enzyme_name;
					}
					
					my $x = ( split("\t", $enzymes_db{$enzyme_name}) )[4];
				}
			}
			
			@selected_enz_names = ();
			foreach my $filtered_enzyme (@filtered_enzymes)
			{
				push @selected_enz_names, $filtered_enzyme;
			}
			
			$R_frame_selEnzymes_listBox->delete(0, 'end');
			foreach my $enzyme_name (@selected_enz_names)
			{
				$R_frame_selEnzymes_listBox->insert('end', $enzyme_name);
			}
			
			$R_frame_selEnzymes_No->configure(-text => scalar(@selected_enz_names) . " enzymes" );
		}
	}
)->pack(-side => 'top', -anchor => 'w', -padx => 5, -pady => 5);

$R_frame_enzymeProperties_label_text = $R_frame->Text(-insertwidth => 0, -cursor => 'left_ptr', -insertontime => 0, -insertofftime => 0, -width => 17, -background => 'gray95', -relief => 'groove', -pady => 5, -padx => 5, -spacing3 => 10, -height => 6);
$R_frame_enzymeProperties_label_text->insert('end', "Enzyme\n");
$R_frame_enzymeProperties_label_text->insert('end', "Recogn. seq.\n");
$R_frame_enzymeProperties_label_text->insert('end', "Cutting site\n");
$R_frame_enzymeProperties_label_text->insert('end', "Overhang\n");
$R_frame_enzymeProperties_label_text->tagConfigure('iso', -foreground => 'black');
$R_frame_enzymeProperties_label_text->tagConfigure('comp', -foreground => 'black');
$R_frame_enzymeProperties_label_text->insert('end', "Isoschizomers\n",'iso');

$R_frame_enzymeProperties_label_text->insert('end', "Companies",'comp');

$R_frame_enzymeProperties_label_text->tagBind('iso', '<Button-1>', sub {
	if (!defined $iso_state or $iso_state == 0) { $iso_state = 1 }
	else { $iso_state = 0 }
	
	$R_frame_enzymeProperties_values_text->delete('0.0', 'end');
	
	if ($iso_state == 1)
	{
		if (( split("\t", $enzymes_db{$selected_enz_name}) )[3] eq "null")
		{
			$R_frame_enzymeProperties_values_text->insert('end', "None");	
		}
		else
		{
			my $first = 0;
			foreach my $isoschisomer ( split(",", ( split("\t", $enzymes_db{$selected_enz_name}) )[3] ) )
			{
				if ($first == 0) { $R_frame_enzymeProperties_values_text->insert('end', "$isoschisomer") }
				else { $R_frame_enzymeProperties_values_text->insert('end', "\n$isoschisomer") }
				$first++;
			}
		}
		$R_frame_enzymeProperties_label_text->tagConfigure('iso', -foreground => 'red');
		$R_frame_enzymeProperties_label_text->tagConfigure('comp', -foreground => 'black');
		$comp_state = 0;
	}
	else
	{
		my @data = split("\t", $enzymes_db{$selected_enz_name});
		my @seq;
		for ( my $i = 0; $i < scalar( split("", $data[1]) ); $i++ )
		{
			if ($data[0] == $i) { push @seq, "'" }
			else
			{
				push @seq, ( split("", $data[1]) )[$i]
			}
		}
		$R_frame_enzymeProperties_values_text->delete('0.0', 'end');
		$R_frame_enzymeProperties_values_text->insert('end', $selected_enz_name . "\n");
		$R_frame_enzymeProperties_values_text->insert('end', $data[-1] . "\n");
		$R_frame_enzymeProperties_values_text->insert('end', $data[0] . "\n");
		$R_frame_enzymeProperties_values_text->insert('end', $data[2]);		
		$R_frame_enzymeProperties_label_text->tagConfigure('iso', -foreground => 'black');
	}
});

$R_frame_enzymeProperties_label_text->tagBind('comp', '<Button-1>', sub {
	if (!defined $comp_state or $comp_state == 0) { $comp_state = 1 }
	else { $comp_state = 0 }
	
	$R_frame_enzymeProperties_values_text->delete('0.0', 'end');
	
	if ($comp_state == 1)
	{
		if (( split("\t", $enzymes_db{$selected_enz_name}) )[4] eq "null")
		{
			$R_frame_enzymeProperties_values_text->insert('end', "None");			
		}
		else
		{
			my $first = 0;
			foreach my $company_ID ( split("", ( split("\t", $enzymes_db{$selected_enz_name}) )[4] ) )
			{
				F: foreach my $company ( split("\t", $enzymes_db{companies}) )
				{
					if ( ( split(",", $company) )[0] eq  $company_ID)
					{
						my $company_name = ( split(",", $company) )[1];
						if ($first == 0) { $R_frame_enzymeProperties_values_text->insert('end', "$company_name") }
						else { $R_frame_enzymeProperties_values_text->insert('end', "\n$company_name") }
						$first++;
						last F;
					}
				}
			}
		}
		$R_frame_enzymeProperties_label_text->tagConfigure('comp', -foreground => 'red');
		$R_frame_enzymeProperties_label_text->tagConfigure('iso', -foreground => 'black');
		$iso_state = 0;
	}
	else
	{
		my @data = split("\t", $enzymes_db{$selected_enz_name});
		my @seq;
		for ( my $i = 0; $i < scalar( split("", $data[1]) ); $i++ )
		{
			if ($data[0] == $i) { push @seq, "'" }
			else
			{
				push @seq, ( split("", $data[1]) )[$i]
			}
		}
		$R_frame_enzymeProperties_values_text->delete('0.0', 'end');
		$R_frame_enzymeProperties_values_text->insert('end', $selected_enz_name . "\n");
		$R_frame_enzymeProperties_values_text->insert('end', $data[-1] . "\n");
		$R_frame_enzymeProperties_values_text->insert('end', $data[0] . "\n");
		$R_frame_enzymeProperties_values_text->insert('end', $data[2]);
		
		$R_frame_enzymeProperties_label_text->tagConfigure('comp', -foreground => 'black');
	}
});
$R_frame_enzymeProperties_values_text = $R_frame->Scrolled('Text', -scrollbars => 'e', -insertwidth => 0, -insertontime => 0, -insertofftime => 0, -width => 17, -wrap => 'word', -background => 'gray95', -relief => 'groove', -pady => 5, -padx => 5, -spacing3 => 10, -height => 6);

my $L_upper_frame = $L_frame->Frame->pack(
	-side => 'top',
	-fill => 'x',
	-padx => 5
);
my $L_center_frame = $L_frame->Frame->pack(
	-side => 'top'
);
my $L_lower_frame = $L_frame->Frame->pack(
	-side => 'top',
	-fill => 'x',
	-padx => 5
);

# $L_upper_1_frame - where '1' means the row number
my $L_upper_1_frame = $L_upper_frame->Frame->pack(-side => 'top', -fill => 'x');
my $L_upper_2_frame = $L_upper_frame->Frame->pack(-side => 'top', -fill => 'x');
my $L_upper_3_frame = $L_upper_frame->Frame->pack(-side => 'top', -fill => 'x');
my $L_upper_4_frame = $L_upper_frame->Frame->pack(-side => 'top', -fill => 'x');
my $L_upper_5_frame = $L_upper_frame->Frame->pack(-side => 'top', -fill => 'x');

# $L_upper_1_2_frame - where '2' means the column number and '1' - the row number
my $L_upper_1_1_frame = $L_upper_1_frame->Frame->pack(-side => 'left');
my $L_upper_1_2_frame = $L_upper_1_frame->Frame->pack(-side => 'left');
my $L_upper_1_3_frame = $L_upper_1_frame->Frame->pack(-side => 'left');
my $L_upper_1_4_frame = $L_upper_1_frame->Frame->pack(-side => 'left');

my $L_upper_2_1_frame = $L_upper_2_frame->Frame->pack(-side => 'left');
my $L_upper_2_2_frame = $L_upper_2_frame->Frame->pack(-side => 'left');
my $L_upper_2_3_frame = $L_upper_2_frame->Frame->pack(-side => 'left');
my $L_upper_2_4_frame = $L_upper_2_frame->Frame->pack(-side => 'left');

my $L_upper_3_1_frame = $L_upper_3_frame->Frame->pack(-side => 'left');
my $L_upper_3_2_frame = $L_upper_3_frame->Frame->pack(-side => 'left');
my $L_upper_3_3_frame = $L_upper_3_frame->Frame->pack(-side => 'left');
my $L_upper_3_4_frame = $L_upper_3_frame->Frame->pack(-side => 'left');

my $L_upper_4_1_frame = $L_upper_4_frame->Frame->pack(-side => 'left');
my $L_upper_4_2_frame = $L_upper_4_frame->Frame->pack(-side => 'left');
my $L_upper_4_3_frame = $L_upper_4_frame->Frame->pack(-side => 'left');
my $L_upper_4_4_frame = $L_upper_4_frame->Frame->pack(-side => 'left');

my $L_upper_5_1_frame = $L_upper_5_frame->Frame->pack(-side => 'left');
my $L_upper_5_2_frame = $L_upper_5_frame->Frame->pack(-side => 'left');
my $L_upper_5_3_frame = $L_upper_5_frame->Frame->pack(-side => 'left');
my $L_upper_5_4_frame = $L_upper_5_frame->Frame->pack(-side => 'left');

$L_upper_1_1_frame->Label(-text => 'Input files', -width => 12)->pack(-side => 'top');
$L_upper_2_1_frame->Label(-text => 'Enzymes', -width => 12)->pack(-side => 'top');
$L_upper_3_1_frame->Label(-text => 'Reference', -width => 12)->pack(-side => 'top');
$L_upper_4_1_frame->Label(-text => 'VCF', -width => 12)->pack(-side => 'top');

$L_upper_1_2_frame->Label(-text => 'Selected file', -width => 23)->pack(-side => 'top');
my $enzyme_entry = $L_upper_2_2_frame->Entry(-insertwidth => 1, -width => 20,-textvariable => \$enzyme_file_name)->pack(-side => 'top',-pady => 1.3);
my $reference_entry = $L_upper_3_2_frame->Entry(-insertwidth => 1, -width => 20,-textvariable => \$reference_file_name)->pack(-side => 'top',-pady => 1.3);
my $raw_vcf_entry = $L_upper_4_2_frame->Entry(-insertwidth => 1, -width => 20,-textvariable => \$raw_vcf_file_name)->pack(-side => 'top',-pady => 1.3);

$L_upper_1_3_frame->Label(-width => 7)->pack(-side => 'top');
my $enzyme_chooseFile_button;
$enzyme_chooseFile_button = $L_upper_2_3_frame->Button(
	-text => 'D',
	-image => $folder_image,
	-command => sub { fileDialog($mw,"enzyme") }
)->pack(
	-side => 'left'
);
my $enzyme_analyze_button;
my $L_lower_col1_mining_button;
$enzyme_analyze_button = $L_upper_2_3_frame->Button(
	-text => 'D',
	-image => $analyze_image,
	-command => sub {
		if (defined $enzyme_file_name and -e $enzyme_file_name and $jobID == 0)
		{
			start_enzymes_check();
			$enzyme_analyze_button->configure(-state => 'disabled');
			$enzyme_chooseFile_button->configure(-state => 'disabled');
		}
		elsif (defined $enzyme_file_name and !-e $enzyme_file_name)
		{
			$enzyme_check->configure(-image => $fail_image);
			curr_time();
			$terminal->insert('end', "Warning", 'warning');
			$terminal->insert('end', " - the file does not exist.\n\n");
			$terminal->see('end');
		}		
	}
)->pack(
	-side => 'left'
);
my $reference_chooseFile_button;
$reference_chooseFile_button = $L_upper_3_3_frame->Button(
	-text => 'D',
	-image => $folder_image,
	-command => sub { fileDialog($mw,"reference") }
)->pack(
	-side => 'left'
);
my $reference_analyze_button;
$reference_analyze_button = $L_upper_3_3_frame->Button(
	-text => 'D',
	-image => $analyze_image,
	-command => sub {
		if (defined $reference_file_name and -e $reference_file_name and $jobID == 0)
		{
			start_reference_check();
			$reference_analyze_button->configure(-state => 'disabled');
			$reference_chooseFile_button->configure(-state => 'disabled');
		}
		elsif (defined $reference_file_name and !-e $reference_file_name)
		{
			$reference_check->configure(-image => $fail_image);
			curr_time();
			$terminal->insert('end', "Warning", 'warning');
			$terminal->insert('end', " - the file does not exist.\n");
			$terminal->see('end');
		}
		

	}
)->pack(
	-side => 'left'
);
my $raw_chooseFile_button;
$raw_chooseFile_button = $L_upper_4_3_frame->Button(
	-text => 'D',
	-image => $folder_image,
	-command => sub { fileDialog($mw,"raw_vcf") }
)->pack(
	-side => 'left'
);
my $raw_vcf_analyze_button;
$raw_vcf_analyze_button = $L_upper_4_3_frame->Button(
	-image => $analyze_image,
	-command => sub {
		if (defined $raw_vcf_file_name and -e $raw_vcf_file_name and $jobID == 0 and $reference_analysis_results[0] == 1)
		{
			raw_start_vcf_check();
			$raw_vcf_analyze_button->configure(-state => 'disabled');
			$raw_chooseFile_button->configure(-state => 'disabled');
		}
		elsif (defined $raw_vcf_file_name and -e $raw_vcf_file_name and $jobID > 0)
		{
			curr_time();
			$terminal->insert('end', "Warning", 'warning');
			$terminal->insert('end', " - please, wait until the current running process is finished.\n");
			$terminal->see('end');
		}
		elsif (defined $raw_vcf_file_name and -e $raw_vcf_file_name and !$reference_analysis_results[0] == 1)
		{
			curr_time();
			$terminal->insert('end', "Warning", 'warning');
			$terminal->insert('end', " - before parsing VCF file, the reference file must first be loaded.\n\n");
			$terminal->see('end');
		}
		elsif (defined $raw_vcf_file_name and !-e $raw_vcf_file_name)
		{
			curr_time();
			$raw_vcf_check->configure(-image => $fail_image);
			$terminal->insert('end', "Warning", 'warning');
			$terminal->insert('end', " - the file does not exist.\n\n");
			$terminal->see('end');
		}		
	}
)->pack(
	-side => 'left'
);

$L_upper_1_4_frame->Label(
	-text => 'Status',
)->pack(
	-side => 'top'
);

$enzyme_check = $L_upper_2_4_frame->Label(
	-text => 'none',
	-foreground => 'grey',
)->pack(
	-side => 'left',
	-padx => 5
);

$enzyme_check_status = $L_upper_2_4_frame->Label(
	-padx => 5
)->pack(
	-side => 'left'
);

$reference_check = $L_upper_3_4_frame->Label(
	-text => 'none',
	-foreground => 'grey'	
)->pack(
	-side => 'left',
	-padx => 5
);

$reference_check_status = $L_upper_3_4_frame->Label(
	-padx => 5
)->pack(
	-side => 'left'
);

$raw_vcf_check = $L_upper_4_4_frame->Label(
	-text => 'none',
	-foreground => 'grey',
)->pack(
	-side => 'left',
	-padx => 5
);

$raw_vcf_check_status = $L_upper_4_4_frame->Label(
	-padx => 5
)->pack(
	-side => 'left'
);


my $L_center_col1_frame = $L_center_frame->Frame->pack(-side => 'left', -fill => 'x', -pady => 10, -padx => 5);
my $L_center_col2_frame = $L_center_frame->Frame->pack(-side => 'left', -fill => 'x');
my $L_center_col1_2_label = $L_center_col1_frame->Label(-anchor => 'w', -text => "DNA sequence length flanking the variant in the output file:")->pack(-side => 'top', -fill => 'x');
$L_center_col2_2_entry = $L_center_col2_frame->Entry(-insertwidth => 1, -width => 5, -textvariable => \$output_seq_len, -justify => 'right')->pack(-side => 'left', -padx => 0);
$L_center_col2_frame->Label(-text => 'bp')->pack(-side => 'top', -padx => 0);
my $L_lower_container_frame = $L_lower_frame->Frame->pack(-side => 'top', -fill => 'x');
my $L_lower_row0_frame = $L_lower_container_frame->Frame->pack(-side => 'top', -fill => 'x');
my $L_lower_row1_frame = $L_lower_container_frame->Frame->pack(-side => 'top', -fill => 'x');
my $L_lower_row3_frame = $L_lower_container_frame->Frame->pack(-side => 'top', -fill => 'x');

$L_lower_col1_mining_button = $L_lower_row1_frame->Button(
	-text => 'Start CAPS mining',
	-width => 15,
	-state => 'disabled',
	-command => sub {
		if ($output_seq_len !~ /^[0-9]+$/ or ($custom == 1 and $custom_value !~ /^[0-9,]+$/) )
		{
			curr_time();
			$terminal->insert('end', "Warning", 'warning');
			$terminal->insert('end', " - something is wrong with the marked parameter value. The parameter must be numerical.\n\n");
			$terminal->see('end');
			if ($output_seq_len !~ /^[0-9]+$/) { $L_center_col2_2_entry->configure(-background => 'red') } else { $L_center_col2_2_entry->configure(-background => 'white') }
			if ($custom == 1 and $custom_value !~ /^[0-9]+$/) { $R_frame_parameters_custom_entry->configure(-background => 'red') } else { $R_frame_parameters_custom_entry->configure(-background => 'white') }
		}
		else
		{
			$L_center_col2_2_entry->configure(-background => 'white');
			$R_frame_parameters_custom_entry->configure(-background => 'white');
			
			if ( scalar(@selected_enz_names) > 0 and $enzyme_analysis_results[0] == 1 and $reference_analysis_results[0] == 1 and $vcf_analysis_results{err_code} == 1 and $jobID == 0)
			{
				$L_lower_col1_mining_button->configure(-state => 'disabled');
				$actualSNPNo = 0;
				start_caps_mining();
			}
			elsif ( scalar(@selected_enz_names) == 0 )
			{
				curr_time();
				$terminal->insert('end', "Warning", 'warning');
				$terminal->insert('end', " - cannot proceed. Please, select some restriction enzymes.\n\n");
				$terminal->see('end');
			}
		}
	}
)->pack(-side => 'left', -anchor => 'w');

my $caps_mining_progress_frame = $L_lower_row1_frame->Frame;
my $caps_mining_result_label = $L_lower_row1_frame->Label->pack(-side => 'left', -anchor => 'w');
my $caps_mining_prepare_enzymes_label = $L_lower_row1_frame->Label();
my $caps_mining_stop_button = $caps_mining_progress_frame->Button(-image => $cancel_image, -command => sub { $stop = 1 } )->pack(-side => 'left', -anchor => 'w');;
my $L_lower_col1_mining_textFrame = $caps_mining_progress_frame->Text(-width => 13, -height => 1, -state => 'disabled')->pack(-side => 'left', -anchor => 'w', -padx => 5);
my $progressBar = $L_lower_col1_mining_textFrame->ProgressBar(-variable => \$capsMining_percent, -width => 14, -length => 90, -gap => 0, -from => 0, -to => 100, -foreground => 'blue', -troughcolor => 'white');
$L_lower_col1_mining_textFrame->windowCreate('end', -window => $progressBar);
my $caps_mining_progress_label = $caps_mining_progress_frame->Label()->pack(-side => 'left', -anchor => 'w');

my $cfw_open_button = $L_lower_row3_frame->Button(
	-text => 'Filtration utilities',
	-width => 15,
	-state => 'normal',
	-command => sub {
		$Caps_Filtration_Window->deiconify;
		$Caps_Filtration_Window->raise;
	}
	)->pack(-side => 'left', -anchor => 'w');

$terminal->pack(-padx => 5, -pady => 5, -expand => 1, -fill => 'both');
$terminal->insert('end', "VCF2CAPS v2.0\n\n");
$terminal->insert('end', "Welcome ...\n\n");

$mw->deiconify;
$mw->raise;

MainLoop;


#-------------------------------------------------#
# The subroutine to chose a new working directory #
#-------------------------------------------------#
sub new_working_directory
{
	my $working_dir_ref = $_[0];
	my $new_working_dir = $mw->chooseDirectory(-initialdir => '.', -title => 'Choose a working directory');
	if ($^O eq 'MSWin32')
	{
		$new_working_dir = Encode::encode("windows-1252", $new_working_dir);
	}

	if (!defined $new_working_dir)
	{
		curr_time();
		$terminal->insert('end', "No directory selected.\n\n");
		$terminal->see('end');
	}
	else
	{
		curr_time();
		$terminal->insert('end', "Selected '");		
		$terminal->insert('end', "$new_working_dir", 'mark');
		$terminal->insert('end', "' as a working directory.\n\n");
		$terminal->see('end');
		$$working_dir_ref = $new_working_dir . "/";
	}
}

#-------------------------------------------------------------#
# The subroutine to save a chosen enzymes list to the text file #
#-------------------------------------------------------------#
sub fileDialog_save_enzyme
{
	my $w = shift;
	my $file;
	my @types =
		(["Text files", '.txt','TEXT'],
		["All files", '*']
	);
	
	$file = $w->getSaveFile(-filetypes => \@types, -initialfile => 'selected_enzymes', -defaultextension => '.txt');
	if (defined $file and $file ne "")
	{
		if (-e $file) { unlink $file }
		if ( scalar(@selected_enz_names) > 0 )
		{
			my $err = 0;
			open my $fh, '>>', $file or $err = 1;
				foreach my $enzyme (@selected_enz_names)
				{
					print $fh "$enzyme\n";
				}
			close $fh;
			
			if ($err == 0)
			{
				curr_time();
				$terminal->insert('end', "The selected enzymes were saved sucessfully to the file '");
				$file =~ s/.*[\\\/]//g;
				$terminal->insert('end', "$file", 'mark');
				$terminal->insert('end', "'.\n\n");
				$terminal->see('end');
			}
			else
			{
				curr_time();
				$terminal->insert('end', "Warning", 'warning');
				$terminal->insert('end', " - there was a problem during saving the file.\n\n");
				$terminal->see('end');
			}
		}
		else
		{
			curr_time();
			$terminal->insert('end', "Warning", 'warning');
			$terminal->insert('end', " - the selected enzymes frame is empty. Please, select some enzymes and then try to save them.\n\n");
			$terminal->see('end');
		}
	}
}

#----------------------------------------------------------------------------#
# The subroutine to retrieve a path to the file containing custom enzyme names #
#----------------------------------------------------------------------------#
sub fileDialog_load_enzyme {
	my $w = shift;
	my $file;
	my @types =
		(["Text files", '.txt','TEXT'],
		["All files", '*']
	);
	
	if ( scalar(@allEnzymesNames) > 0 )
	{	
		$file = $w->getOpenFile(-filetypes => \@types);
		if (defined $file and $file ne "")
		{
			my @loaded_enzymes;
			my $err = 0;
			my $fh;
			
			if ( !open $fh, '<', $file )
			{
				my $err_code = $!;
				my $file_tmp = $file;
				$file_tmp =~ s/.*[\\\/]//g;
				curr_time();
				$terminal->insert('end', "Warning", 'warning');
				$terminal->insert('end', " - there was a problem during opening the file: '");
				$terminal->insert('end', "$file_tmp", 'mark');
				$terminal->insert('end', "'.\n\n");
				$terminal->see('end');				
				warning_dial("Cannot open the file ::: $err_code");
				return;
			}
			
			while (<$fh>)
			{
				chomp $_;
				if ($_ =~ /^[A-Za-z0-9-]+$/) { push @loaded_enzymes, $_ }
				else { $err = 2; last }
			}
			
			close $fh;
			
			if ($err == 2)
			{
				curr_time();
				$terminal->insert('end', "Warning", 'warning');
				$terminal->insert('end', " - something is wrong with the file content. Does it really contain the enzymes list?\n\n");
				$terminal->see('end');				
			}
			else
			{
				my @loaded_enzymes_OK;
				my @loaded_enzymes_fail = @loaded_enzymes;
				foreach my $loaded_enzyme (@loaded_enzymes)
				{
					foreach my $enzyme (@allEnzymesNames)
					{
						if ($enzyme eq $loaded_enzyme)
						{
							push @loaded_enzymes_OK, $loaded_enzyme;
							@loaded_enzymes_fail = grep { $_ ne $loaded_enzyme } @loaded_enzymes_fail;
							last;
						}
					}
				}
				
				@selected_enz_names = ();
				foreach my $OK (@loaded_enzymes_OK)
				{
					push @selected_enz_names, $OK;
				}
				
				$R_frame_selEnzymes_listBox->delete(0, 'end');
				foreach my $enzyme_name (@selected_enz_names)
				{
					$R_frame_selEnzymes_listBox->insert('end', $enzyme_name);
				}
				$R_frame_selEnzymes_No->configure(-text => scalar(@selected_enz_names) . " enzymes" );
				
				if ( scalar(@loaded_enzymes_fail) > 0 )
				{
					my $fail_No = scalar(@loaded_enzymes_fail);
					curr_time();
					$terminal->insert('end', "Warning", 'warning');
					$terminal->insert('end', " - $fail_No", 'mark');
					if ($fail_No == 1)
					{
						$terminal->insert('end', " enzyme is not present in the database: ");
						$terminal->insert('end', "$loaded_enzymes_fail[0]", 'mark');
						$terminal->insert('end', ". The enzyme was discarded from the list.\n\n");
						$terminal->see('end');						
					}
					else
					{
						$terminal->insert('end', " enzymes are not present in the database: ");
						for (my $i = 0; $i < scalar(@loaded_enzymes_fail); $i++)
						{
							if ($i == scalar(@loaded_enzymes_fail) - 1 )
							{
								$terminal->insert('end', "$loaded_enzymes_fail[$i]", 'mark');
							}
							else
							{
								$terminal->insert('end', "$loaded_enzymes_fail[$i]", 'mark');
								$terminal->insert('end', ", ");
							}
						}
						$terminal->insert('end', ". The enzymes were discarded from the list.\n\n");
						$terminal->see('end');
					}
				}
				
				if ( scalar(@loaded_enzymes_OK) > 0 )
				{
					my $OK_No = scalar(@loaded_enzymes_OK);
					curr_time();
					$terminal->insert('end', "$OK_No", 'mark');
					if ( scalar(@loaded_enzymes_OK) == 1 )
					{
						$terminal->insert('end', " selected enzyme was loaded sucessfully.\n\n");
						$terminal->see('end');
					}
					else
					{
						$terminal->insert('end', " selected enzymes were loaded sucessfully.\n\n");
						$terminal->see('end');
					}
				}
			}
		}
		else
		{
			warning_dial('Warning ::: The file could not be loaded.');
		}		
	}
	else
	{
		curr_time();
		$terminal->insert('end', "Warning", 'warning');
		$terminal->insert('end', " - please, load the enzyme database before loading selected enzymes.\n\n");
		$terminal->see('end');
	}
}

#-------------------------------------------------------#
# The subroutine opens the URL in the cutom web browser #
#-------------------------------------------------------#
sub open_hyperlink
{
	my $url = shift;
	my $platform = $^O;
	my $cmd;
	if ($platform eq 'darwin') { $cmd = "open \"$url\"" } # Mac OS X
	elsif ($platform eq 'linux') { $cmd = "x-www-browser \"$url\"" } # Linux
	elsif ($platform eq 'MSWin32') { $cmd = "start $url" } # Win95..Win7
	if (defined $cmd) { system($cmd) }
}

#-------------------------------------------------------------#
# The subroutine writes custom messages to the file 'log.txt' #
#-------------------------------------------------------------#
sub LOG
{
	my $text = $_[0];
	
	if ($log_first_use == 0)
	{
		if (-e "log.txt") { unlink "log.txt" }
		
		$log_first_use = 1;
		open my $Ofh, '>>', $working_dir . "log.txt" or warning_dial("Cannot open the file for writing ::: $!");
			my $date = localtime();
			my $date_length = split("", $date);
			print $Ofh "\n#" . "-" x ($date_length + 2) . "#\n";
			print $Ofh "# VCF2CAPS v2.0" . " " x (($date_length + 4) - 16) . "#\n";
			print $Ofh "# " . $date . " #\n";
			print $Ofh "#" . "-" x ($date_length + 2) . "#\n";
		close $Ofh;
	}
	
	open my $Ofh, '>>', $working_dir . "log.txt" or warning_dial("Cannot open the file for writing ::: $!");
		print $Ofh $text . "\n";
	close $Ofh;
}

#----------------------------------------------#
# The subroutine downloads the enzyme database #
#----------------------------------------------#
sub start_download_enzyme_db
{
	$download_enzyme_db_result = 0;
	curr_time();
	$terminal->insert('end', "Downloading the database from ");
	$terminal->insert('end', "http://rebase.neb.com/rebase/link_gcg", 'mark');
	$terminal->insert('end', " ...\n\n");
	$terminal->see('end');
	
	my $repeat;
	$jobID = 14;
	$repeat = $mw->repeat( 100 => sub {
		if ( $download_enzyme_db_result != 0 )
		{
			if ( $download_enzyme_db_result == 1 )
			{
				my $enzyme_file_name_tmp = $enzyme_file_name;
				$enzyme_entry->delete(0, 'end');
				$enzyme_entry->insert(0, $enzyme_file_name_tmp);
				$enzyme_entry->xview('end');			
				curr_time();
				$terminal->insert('end', "The database file downloaded sucessfully. Loading the database ...\n\n");
				$terminal->see('end');
				
				if (defined $enzyme_file_name and -f $enzyme_file_name and $jobID == 0)
				{
					start_enzymes_check(); $enzyme_analyze_button->configure(-state => 'disabled');
				}
				elsif (defined $enzyme_file_name and !-f $enzyme_file_name)
				{
					$enzyme_check->configure(-image => $fail_image);
					curr_time();
					$terminal->insert('end', "Warning", 'warning');
					$terminal->insert('end', " - the file does not exist.\n\n");
					$terminal->see('end');
				}
			}
			elsif ( $download_enzyme_db_result == 2 )
			{
				curr_time();
				$terminal->insert('end', "Warning", 'warning');
				$terminal->insert('end', " - something went wrong during download.\n\n");
				$terminal->see('end');
			}
		
			$repeat->cancel;
		}
	} );
}

#--------------------------------------------------------------#
# The subroutine checks the format of the VCF2CAPS output file #
#--------------------------------------------------------------#
sub start_vcf2capsOutput_check
{
	my $filtration_type = shift;
	
	%vcf2capsOutput_results = ();
	$vcf2capsOutput_results{err_code} = 0;
	
	my $repeat;
	if ($filtration_type eq 'scf')
	{
		$jobID = 10;
		$cfw_scf_inputFile_check->configure(-image => $processing_gif);
		$processing_gif->start_animation;
		$repeat = $mw->repeat( 1000 => sub {

			if ( $vcf2capsOutput_results{err_code} != 0 )
			{
				$processing_gif->stop_animation;
				
				if ( $vcf2capsOutput_results{err_code} == 1 )
				{
					$cfw_scf_inputFile_check->configure(-image => $ok_image);
				}
				elsif ( $vcf2capsOutput_results{err_code} == 6 )
				{
					$cfw_scf_inputFile_check->configure(-image => $warning_image);
					my $cfw_scf_input_file_tmp = $cfw_scf_input_file;
					$cfw_scf_input_file_tmp =~ s/.*[\\\/]//g;
					
					curr_time();
					$terminal->insert('end', "Warning", 'warning');
					$terminal->insert('end', " - cannot open the file: '");
					$terminal->insert('end', "$cfw_scf_input_file_tmp", 'mark');
					$terminal->insert('end', "'\n\n");
					$terminal->see('end');
					
					$repeat->cancel;
					warning_dial("Cannot open the file ::: $vcf2capsOutput_results{err_value}");
				}
				else
				{
					$cfw_scf_inputFile_check->configure(-image => '', -text => 'Something is wrong with the file format', -foreground => 'red');
				}
				
				$cfw_scf_inputFile_analyze_button->configure(-state => 'normal');
				$cfw_scf_inputFile_chooseFile_button->configure(-state => 'normal');
				$cfw_scf_start_button->configure(-state => 'normal');
				
				$repeat->cancel;				
			}
		} );
	}
	elsif ($filtration_type eq 'gf')
	{
		$jobID = 9;
		$cfw_gf_inputFile_check->configure(-image => $processing_gif);
		$processing_gif->start_animation;
		$repeat = $mw->repeat( 100 => sub {
			if ( $vcf2capsOutput_results{err_code} != 0 )
			{
				$processing_gif->stop_animation;
				
				if ( $vcf2capsOutput_results{err_code} == 1 )
				{					
					$cfw_gf_inputFile_check->configure(-image => $ok_image);
				}
				elsif ( $vcf2capsOutput_results{err_code} == 6 )
				{
					$cfw_gf_inputFile_check->configure(-image => $warning_image);
					
					my $cfw_gf_input_file_tmp = $cfw_gf_input_file;
					$cfw_gf_input_file_tmp =~ s/.*[\\\/]//g;
					
					curr_time();
					$terminal->insert('end', "Warning", 'warning');
					$terminal->insert('end', " - cannot open the file: '");
					$terminal->insert('end', "$cfw_gf_input_file_tmp", 'mark');
					$terminal->insert('end', "'\n\n");
					$terminal->see('end');
					
					$repeat->cancel;
					
					warning_dial("Cannot open the file ::: $vcf2capsOutput_results{err_value}");
				}
				else
				{
					$cfw_gf_inputFile_check->configure(-image => '', -text => 'Something is wrong with the file format', -foreground => 'red');
				}
				
				$cfw_gf_inputFile_analyze_button->configure(-state => 'normal');
				$cfw_gf_inputFile_chooseFile_button->configure(-state => 'normal');
				$cfw_gf_start_button->configure(-state => 'normal');
				
				$repeat->cancel;
				
			}
		} );
	}
	elsif ($filtration_type eq 'c2f')
	{
		$jobID = 11;
		$cfw_c2f_inputFile_check->configure(-image => $processing_gif);
		$processing_gif->start_animation;
		$repeat = $mw->repeat( 100 => sub {
			if ( $vcf2capsOutput_results{err_code} != 0 )
			{
				$processing_gif->stop_animation;
				
				if ( $vcf2capsOutput_results{err_code} == 1 )
				{					
					$cfw_c2f_inputFile_check->configure(-image => $ok_image);
					
				}
				elsif ( $vcf2capsOutput_results{err_code} == 6 )
				{
					$cfw_c2f_inputFile_check->configure(-image => $warning_image);
					
					my $cfw_c2f_input_file_tmp = $cfw_c2f_input_file;
					$cfw_c2f_input_file_tmp =~ s/.*[\\\/]//g;
					
					curr_time();
					$terminal->insert('end', "Warning", 'warning');
					$terminal->insert('end', " - cannot open the file: '");
					$terminal->insert('end', "$cfw_c2f_input_file_tmp", 'mark');
					$terminal->insert('end', "'\n\n");
					$terminal->see('end');
					
					$repeat->cancel;
					warning_dial("Cannot open the file ::: $vcf2capsOutput_results{err_value}");
					
				}
				else
				{
					$cfw_c2f_inputFile_check->configure(-image => '', -text => 'Something is wrong with the file format', -foreground => 'red');
				}
				
				$cfw_c2f_inputFile_analyze_button->configure(-state => 'normal');
				$cfw_c2f_inputFile_chooseFile_button->configure(-state => 'normal');
				$cfw_c2f_start_button->configure(-state => 'normal');
				
				$repeat->cancel;				
			}
		} );
	}	
}

#---------------------------------------------------------------------------------------#
# The subroutine triggers and checks the progress of enzyme database file checking step #
#---------------------------------------------------------------------------------------#
sub start_enzymes_check
{
	if (defined $enzyme_file_name and -e $enzyme_file_name)
	{
		@enzyme_analysis_results = (0);
		$jobID = 2;
		my $repeat;
		my $c = 0;
		$enzyme_check_status->configure(-text => "");
		$enzyme_check->configure(-image => $processing_gif);
		$processing_gif->start_animation;
		$repeat = $mw->repeat( 100 => sub {
			if ($enzyme_analysis_results[0] != 0)
			{
				$processing_gif->stop_animation;
				
				if ($enzyme_analysis_results[0] == 1)
				{
					$enzyme_check->configure(-image => $ok_image);
					
					$R_frame_allEnzymes_listBox->delete(0,'end');
					
					foreach my $enzyme_name (@allEnzymesNames)
					{
						$R_frame_allEnzymes_listBox->insert('end', $enzyme_name);
					}
					
					$R_frame_allEnzymes_No->configure(-text => scalar(@allEnzymesNames) . " enzymes" );
					curr_time();
					$terminal->insert('end', "Loaded REBASE database (version ");
					$terminal->insert('end', $enzymes_db{date}, 'mark');
					$terminal->insert('end', ") comprising ");
					$terminal->insert('end', scalar(@allEnzymesNames), 'mark');
					$terminal->insert('end', " enzymes.\n\n");
					$terminal->see('end');
					
					if ($reference_analysis_results[0] == 1 and $enzyme_analysis_results[0] == 1 and $vcf_analysis_results{err_code} == 1)
					{
						$L_lower_col1_mining_button->configure(-state => 'normal');
					}
				}
				elsif ($enzyme_analysis_results[0] == 3)
				{
					$enzyme_check->configure(-image => $fail_image);

					curr_time();
					$terminal->insert('end', "Warning", 'warning');
					$terminal->insert('end', " - something is wrong with the database file.\n\n");
					$terminal->see('end');
				}
				elsif ($enzyme_analysis_results[0] == 4)
				{
					$enzyme_check->configure(-image => $fail_image);
					
					curr_time();
					$terminal->insert('end', "Warning", 'warning');
					$terminal->insert('end', " - cannot open the enzyme database file.\n\n");
					$terminal->see('end');	
				}
				
				$enzyme_analyze_button->configure(-state => 'normal');
				$enzyme_chooseFile_button->configure(-state => 'normal');
				$repeat->cancel;
				
				if ( $enzyme_analysis_results[0] == 4 )
				{
					warning_dial("Cannot open the file ::: $enzyme_analysis_results[1]");
				}
			}
		} );
	}
	elsif (defined $enzyme_file_name)
	{
		$enzyme_check->configure(-image => $fail_image);

		curr_time();
		$terminal->insert('end', "Warning", 'warning');
		$terminal->insert('end', " - the file does not exist.\n\n");
		$terminal->see('end');
	}	
}

#---------------------------------------------------------------------------------#
# The subroutine triggers and checks the progress of reference file checking step #
#---------------------------------------------------------------------------------#
sub start_reference_check
{
	if (defined $reference_file_name and -e $reference_file_name)
	{
		my $reference_file_name_tmp;
		$reference_file_name_tmp = $reference_file_name;
		$reference_file_name_tmp =~ s/.*[\\\/]//g;	
		
		my $reference_index_file_name = $reference_file_name;
		$reference_index_file_name =~ s/\..+$/\.index/g;
		$reference_index_file_name =~ s/.*[\\\/]//g;				
		
		%reference_index_check_result = ();
		$reference_index_check_result{error_code} = "";
		$jobID = 13;
		
		curr_time();
		$terminal->insert('end', "Checking MD5 of the reference file, please wait ...\n\n",);
		$terminal->see('end');
		
		$reference_check->configure(-image => $processing_gif);
		$processing_gif->start_animation;
		
		my $repeat_reference_check;
		$repeat_reference_check = $mw->repeat( 100 => sub {
			if ( $reference_index_check_result{error_code} ne "" )
			{
				if ( $reference_index_check_result{error_code} eq "OK" )
				{
					$processing_gif->stop_animation;
					$reference_check->configure(-image => $ok_image);
					curr_time();
					$terminal->insert('end', "Using previously created index file: '");
					$terminal->insert('end', $reference_index_file_name, 'mark',);
					$terminal->insert('end', "'\n\n");
					$terminal->see('end');
					
					if ($seqExtractor_error_code[1] == 1 and $reference_analysis_results[0] == 1 and $enzyme_analysis_results[0] == 1 and $vcf_analysis_results{err_code} == 1)
					{
						$L_lower_col1_mining_button->configure(-state => 'normal');
					}
				}
				else
				{
					if ( $reference_index_check_result{error_code} eq "wrong_md5" )
					{
						curr_time();
						$terminal->insert('end', "The index file:  '");
						$terminal->insert('end', $reference_index_file_name, 'mark',);
						$terminal->insert('end', "' in the working direcotry, is related to the different reference file. \n");
						$terminal->insert('end', "Deleting '");
						$terminal->insert('end', $reference_index_file_name, 'mark',);
						$terminal->insert('end', "' ...\n\n");
						$terminal->see('end');
						
						unlink $working_dir . $reference_index_file_name;
					}
					elsif ( $reference_index_check_result{error_code} eq "no_md5" )
					{
						curr_time();
						$terminal->insert('end', "The index file:  '");
						$terminal->insert('end', $reference_index_file_name, 'mark',);
						$terminal->insert('end', "' in the working direcotry, does not have the md5 hash header.\n");
						$terminal->insert('end', "Deleting '");
						$terminal->insert('end', $reference_index_file_name, 'mark',);
						$terminal->insert('end', "' ...\n\n");
						$terminal->see('end');
					}
					elsif ( $reference_index_check_result{error_code} eq "cannot_open_ref" )
					{
						$processing_gif->stop_animation;
						$reference_check->configure(-image => $warning_image);
						
						curr_time();
						$terminal->insert('end', "Warning", 'warning');
						$terminal->insert('end', " - cannot open the file: '");
						$terminal->insert('end', $reference_file_name_tmp, 'mark',);
						$terminal->insert('end', "'\n\n");
						$terminal->see('end');
						
						$L_lower_col1_mining_button->configure(-state => 'disabled');
						
						$repeat_reference_check->cancel;
						
						warning_dial("Cannot open the file ::: $reference_index_check_result{error_value}");
						return;
					}
					elsif ( $reference_index_check_result{error_code} eq "cannot_open_index" )
					{
						$processing_gif->stop_animation;
						$reference_check->configure(-image => $warning_image);
						
						curr_time();
						$terminal->insert('end', "Warning", 'warning');
						$terminal->insert('end', " - cannot open the file: '");
						$terminal->insert('end', $reference_index_file_name, 'mark',);
						$terminal->insert('end', "'\n\n");
						$terminal->see('end');
						
						$L_lower_col1_mining_button->configure(-state => 'disabled');
						
						$repeat_reference_check->cancel;
						
						warning_dial("Cannot open the file ::: $reference_index_check_result{error_value}");
						return;
					}
					
					if ( -e $working_dir . $reference_index_file_name) { unlink $working_dir . $reference_index_file_name }
					
					curr_time();
					$terminal->insert('end', "Start reference file '");
					$terminal->insert('end', "$reference_file_name_tmp", 'mark');
					$terminal->insert('end', "' integrity check ...\n\n");
					$terminal->see('end');
					
					@reference_analysis_results = (0);
					$jobID = 1;
					
					my $repeat;
					my $c = 0;
					$reference_check->configure(-image => $processing_gif);
					$processing_gif->start_animation;
					$repeat = $mw->repeat( 100 => sub {
						if ($reference_analysis_results[0] != 0)
						{
							$processing_gif->stop_animation;
							
							if ($reference_analysis_results[0] == 1)
							{
								$reference_check->configure(-image => $ok_image);
							
								curr_time();
								$terminal->insert('end', "Integrity of the reference file '");
								$terminal->insert('end', "$reference_file_name_tmp", 'mark');
								$terminal->insert('end', "' confirmed. The index file '");
								$reference_file_name_tmp =~ s/(?<=\.).+$//;
								$terminal->insert('end', "$reference_file_name_tmp" . "index", 'mark');
								$terminal->insert('end', "' was created.\n\n");
								$terminal->see('end');
								
								if ($seqExtractor_error_code[1] == 1 and $reference_analysis_results[0] == 1 and $enzyme_analysis_results[0] == 1 and $vcf_analysis_results{err_code} == 1)
								{
									$L_lower_col1_mining_button->configure(-state => 'normal');
								}
							}
							elsif ($reference_analysis_results[0] == 3)
							{
								$reference_check->configure(-image => $fail_image);

								curr_time();
								$terminal->insert('end', "Warning", 'warning');
								$terminal->insert('end', " - duplicated sequence: ");
								$terminal->insert('end', "$reference_analysis_results[1].\n\n", 'mark');
								$terminal->see('end');
							}
							elsif ($reference_analysis_results[0] == 4)
							{
								$reference_check->configure(-image => $fail_image);
								my @err_data = split(",", $reference_analysis_results[1]);
								curr_time();
								$terminal->insert('end', "Warning", 'warning');
								$terminal->insert('end', " - the refernce file '");
								$terminal->insert('end', "$reference_file_name_tmp", 'mark');
								$terminal->insert('end', "' has different line length in '");
								$terminal->insert('end', "$err_data[0]", 'mark');
								$terminal->insert('end', "' at line ");
								$terminal->insert('end', "$err_data[1]", 'mark');
								$terminal->insert('end', ".\n\n");
								$terminal->see('end');
							}
							elsif ($reference_analysis_results[0] == 5)
							{
								$reference_check->configure(-image => $fail_image);
								
								curr_time();
								$terminal->insert('end', "Warning", 'warning');
								$terminal->insert('end', " - cannot open the file '");
								$terminal->insert('end', "$reference_file_name_tmp", 'mark');
								$terminal->insert('end', "'.\n\n");
								$terminal->see('end');								
							}
							elsif ($reference_analysis_results[0] == 6)
							{
								$reference_check->configure(-image => $fail_image);
								
								curr_time();
								$terminal->insert('end', "Warning", 'warning');
								$terminal->insert('end', " - cannot create the file '");
								$terminal->insert('end', "$reference_file_name_tmp.index", 'mark');
								$terminal->insert('end', "'.\n\n");
								$terminal->see('end');								
							}
							
							$reference_analyze_button->configure(-state => 'normal');
							$reference_chooseFile_button->configure(-state => 'normal');
							$repeat->cancel;
							
							if ( $reference_analysis_results[0] == 5 )
							{
								warning_dial("Cannot open the file ::: $reference_analysis_results[1]");
							}
							elsif ( $reference_analysis_results[0] == 6 )
							{
								warning_dial("Cannot create the file ::: $reference_analysis_results[1]");
							}
						}
					} );
				}
				
				$reference_analyze_button->configure(-state => 'normal');
				$reference_chooseFile_button->configure(-state => 'normal');
				$repeat_reference_check->cancel;
			}
		} );	
	}
	elsif (defined $reference_file_name)
	{
		$reference_check->configure(-image => $fail_image);
		
		curr_time();
		$terminal->insert('end', "Warning", 'warning');
		$terminal->insert('end', " - the file does not exist.\n\n");
		$terminal->see('end');
	}
}


#-------------------------------------------------------------------------------------------------#
# The subroutine triggers and checks the progress of the VCF2CAPS output file to FASTA conversion #
#-------------------------------------------------------------------------------------------------#
sub start_caps_to_fasta_convertion
{
	$cfw_c2f_error_label->packForget;
	$cfw_c2f_convertion_percent = 0;
	@caps_to_fasta_result = (0, 0);
	$caps_filtered = 0;
	
	$cfw_c2f_convertion_percent = ( $caps_filtered / $total_caps_number ) * 100;
	$cfw_c2f_progress_label->configure(-text => sprintf ("%d/%d   %.1f%%", $caps_filtered,$total_caps_number,$cfw_c2f_convertion_percent) );
	$cfw_c2f_progress_frame->pack(-side => 'left', -padx => 5, -pady => 5, -anchor => 'w');
	
	curr_time();
	$terminal->insert('end', "Start convertion of vcf2caps output file into FASTA format ...\n\n");
	$terminal->see('end');
	
	$jobID = 12;
	
	my $repeat;
	$repeat = $mw->repeat( 100 => sub {
		$cfw_c2f_convertion_percent = ( $caps_filtered / $total_caps_number ) * 100;
		
		if ( $caps_to_fasta_result[0] != 0 )
		{	
			my $cfw_c2f_output_file_tmp = $cfw_c2f_output_file;
			$cfw_c2f_output_file_tmp =~ s/.+\/(.+\.fasta)/$1/;
			$cfw_c2f_progress_frame->packForget;
				
			if ( $caps_to_fasta_result[0] == 1 )
			{
				$cfw_c2f_convertion_percent = ( $caps_filtered / $total_caps_number ) * 100;
				$cfw_c2f_progress_label->configure(-text => sprintf ("%d/%d   %.1f%%", $caps_filtered,$total_caps_number,$cfw_c2f_convertion_percent) );
				
				
				
				$cfw_c2f_error_label->configure(-text => "$caps_to_fasta_result[1] sequences were saved to the file: '$cfw_c2f_output_file_tmp'", -foreground => 'black');
				$cfw_c2f_error_label->pack(-side => 'left', -padx => 5, -pady => 5, -anchor => 'w');
				
				curr_time();
				$terminal->insert('end', "Finished convertion of VCF2CAPS output file into FASTA format.");
				$terminal->insert('end', " $caps_to_fasta_result[1]", 'mark');
				$terminal->insert('end', " sequences were saved to the file: '");
				$terminal->insert('end', "$cfw_c2f_output_file_tmp", 'mark');
				$terminal->insert('end', "' file\n\n");
				$terminal->see('end');
			}
			elsif ( $caps_to_fasta_result[0] == 2 )
			{
				$cfw_c2f_progress_frame->packForget;
				
				curr_time();
				$terminal->insert('end', "Warning",'warning');
				$terminal->insert('end', " - cannot write to the output file: ");
				$terminal->insert('end', "$cfw_c2f_output_file_tmp", 'mark');
				$terminal->insert('end', "'.\n\n");
				$terminal->see('end');
			}
			elsif ( $caps_to_fasta_result[0] == 3 )
			{
				$cfw_c2f_progress_frame->packForget;
				
				my $cfw_c2f_input_file_tmp = $cfw_c2f_input_file;
				$cfw_c2f_input_file_tmp =~ s/.*[\\\/]//g;
				$cfw_c2f_inputFile_check->configure(-image => $warning_image);
				
				curr_time();
				$terminal->insert('end', "Warning",'warning');
				$terminal->insert('end', " - cannot open the file: ");
				$terminal->insert('end', "$cfw_c2f_input_file_tmp", 'mark');
				$terminal->insert('end', "'.\n\n");
				$terminal->see('end');
			}
			
			$cfw_c2f_start_button->configure(-state => 'normal');	
			$repeat->cancel;
			
			if ( $caps_to_fasta_result[0] == 2 )
			{
				warning_dial("Cannot write to the file ::: $caps_to_fasta_result[1]");
			}
			elsif ( $caps_to_fasta_result[0] == 3 )
			{
				warning_dial("Cannot open the file ::: $caps_to_fasta_result[1]");
			}
		}		
		elsif ( $caps_filtered > 0 and $total_caps_number > 0 )
		{
			$cfw_c2f_convertion_percent = ( $caps_filtered / $total_caps_number ) * 100;
			$cfw_c2f_progress_label->configure(-text => sprintf ("%d/%d   %.1f%%", $caps_filtered,$total_caps_number,$cfw_c2f_convertion_percent) );
		}
	} );
}


#----------------------------------------------------------------------------------------#
# The subroutine triggers and checks the progress of CAPS markers filtration by genotype #
#----------------------------------------------------------------------------------------#
sub start_caps_filtration
{	
	$cfw_gf_error_label->packForget;
	$cfw_gf_CAPS_filtering_percent = 0;
	$caps_filtration_result[0] = 0;
	$caps_filtration_result[1] = 0;
	$caps_filtered = 0;
	$cfw_gf_start_button->configure(-state => 'disabled');
	my %cfw_groups_tmp;
	my $number_of_filters = 0;
	$cfw_groups{1}{indv} = "";
	$cfw_groups{2}{indv} = "";
	$cfw_groups{3}{indv} = "";
	my $cfw_group_1_value = $cfw_gf_group_1_text->get('1.0', 'end-1c');
	if ($cfw_group_1_value ne "")
	{
		$cfw_groups{1}{max_err} = $cfw_gf_group_1_maxError_value;
		$cfw_groups{1}{indv} = $cfw_group_1_value;
		$number_of_filters++;
	}
	
	my $cfw_group_2_value = $cfw_gf_group_2_text->get('1.0', 'end-1c');
	if ($cfw_group_2_value ne "")
	{
		$cfw_groups{2}{max_err} = $cfw_gf_group_2_maxError_value;
		$cfw_groups{2}{indv} = $cfw_group_2_value;
		$number_of_filters++;
	}
	
	my $cfw_group_3_value = $cfw_gf_group_3_text->get('1.0', 'end-1c');
	if ($cfw_group_3_value ne "")
	{
		$cfw_groups{3}{max_err} = $cfw_gf_group_3_maxError_value;
		$cfw_groups{3}{indv} = $cfw_group_3_value;
		$number_of_filters++;
	}
	
	if ( $number_of_filters <= 1 )
	{
		$cfw_gf_error_label->configure(-text => 'You must input at least two groups', -foreground => 'red');
		$cfw_gf_error_label->pack(-side => 'left', -padx => 5, -pady => 5, -anchor => 'w');
		$cfw_gf_start_button->configure(-state => 'normal');
		return 1;
	}

	$jobID = 8;
	$cfw_gf_CAPS_filtering_percent = ( $caps_filtered / $total_caps_number ) * 100;
	$cfw_gf_progress_label->configure(-text => sprintf ("%d/%d   %.1f%%", $caps_filtered,$total_caps_number,$cfw_gf_CAPS_filtering_percent) );
	$cfw_gf_progress_frame->pack(-side => 'left', -padx => 5, -pady => 5, -anchor => 'w');
	
	curr_time();
	$terminal->insert('end', "Start CAPS markers filtration by genotype ...\n\n");
	$terminal->see('end');
	
	my $repeat;
	$repeat = $mw->repeat( 100 => sub {		
		$cfw_gf_CAPS_filtering_percent = ( $caps_filtered / $total_caps_number ) * 100;
		
		if ( $caps_filtration_result[0] != 0 )
		{
			my $cfw_gf_output_file_tmp = $cfw_gf_output_file;
			$cfw_gf_output_file_tmp =~ s/.+\/(.+\.txt)/$1/;
			
			if ( $caps_filtration_result[0] == 1 )
			{
				$cfw_gf_CAPS_filtering_percent = ( $caps_filtered / $total_caps_number ) * 100;
				$cfw_gf_progress_label->configure(-text => sprintf ("%d/%d   %.1f%%", $caps_filtered,$total_caps_number,$cfw_gf_CAPS_filtering_percent) );
				
				$cfw_gf_progress_frame->packForget;
				
				if ( $caps_filtration_result[1] > 0 )
				{
					$cfw_gf_error_label->configure(-text => "$caps_filtration_result[1] filtered markers were saved to the file: '$cfw_gf_output_file_tmp'", -foreground => 'black');
					$cfw_gf_error_label->pack(-side => 'left', -padx => 5, -pady => 5, -anchor => 'w');
					
					curr_time();
					$terminal->insert('end', "CAPS markers filtration by genotype finished.");
					$terminal->insert('end', " $caps_filtration_result[1]", 'mark');
					$terminal->insert('end', " CAPS were saved to the file: '");
					$terminal->insert('end', "$cfw_gf_output_file_tmp", 'mark');
					$terminal->insert('end', "' file\n\n");
					$terminal->see('end');
				}
				else
				{
					$cfw_gf_error_label->configure(-text => "No CAPS markers were found that meet the specified criteria.", -foreground => 'black');
					$cfw_gf_error_label->pack(-side => 'left', -padx => 5, -pady => 5, -anchor => 'w');
					
					curr_time();
					$terminal->insert('end', "CAPS markers filtration by genotype finished.\n");
					$terminal->insert('end', "Warning",'warning');
					$terminal->insert('end', " - no CAPS markers were found that meet the specified criteria.\n\n");
					$terminal->see('end');
				}
			}
			elsif ( $caps_filtration_result[0] == 2 )
			{
				$cfw_gf_progress_frame->packForget;
				
				curr_time();
				$terminal->insert('end', "Cannot write to the output file: '");				
				$terminal->insert('end', "$cfw_gf_output_file_tmp", 'mark');
				$terminal->insert('end', "' file\n\n");
				$terminal->see('end');
			}
			elsif ( $caps_filtration_result[0] == 3 )
			{
				$cfw_gf_progress_frame->packForget;
				$cfw_gf_inputFile_check->configure(-image => $warning_image);
				
				my $cfw_gf_input_file_tmp = $cfw_gf_input_file;
				$cfw_gf_input_file_tmp =~ s/.*[\\\/]//g;
				
				curr_time();
				$terminal->insert('end', "Cannot open the file: '");				
				$terminal->insert('end', "$cfw_gf_input_file_tmp", 'mark');
				$terminal->insert('end', "' file\n\n");
				$terminal->see('end');
			}
			
			$cfw_gf_start_button->configure(-state => 'normal');						
			$repeat->cancel;
			
			if ( $caps_filtration_result[0] == 2 )
			{
				warning_dial("Cannot write to the file ::: $caps_filtration_result[1]");
			}
			elsif ( $caps_filtration_result[0] == 3 )
			{
				warning_dial("Cannot open the file ::: $caps_filtration_result[1]");
			}
		}		
		elsif ( $caps_filtered > 0 and $total_caps_number > 0 )
		{
			$cfw_gf_CAPS_filtering_percent = ( $caps_filtered / $total_caps_number ) * 100;
			$cfw_gf_progress_label->configure(-text => sprintf ("%d/%d   %.1f%%", $caps_filtered,$total_caps_number,$cfw_gf_CAPS_filtering_percent) );
		}
	} );
}


#----------------------------------------------------------------------------------------#
# The subroutine triggers and checks the progress of the VCF to v2c file convertion step #
#----------------------------------------------------------------------------------------#
sub raw_start_vcf_check
{
	%vcf_analysis_results = ();
	%vcf_analysis_results = (err_code => 0);
	%v2c_check_result = ();
	%v2c_check_result = (err_type => "");
	$jobID = 7;
	my $c = 0;
	$raw_vcf_check->configure(-image => $processing_gif);
	$processing_gif->start_animation;
	
	my $raw_vcf_file_name_tmp = $raw_vcf_file_name;
	$raw_vcf_file_name_tmp =~ s/.*[\\\/]//g;
	my $v2c_file_name_tmp = $raw_vcf_file_name_tmp;
	$v2c_file_name_tmp =~ s/\..+$/\.v2c/;

	my $reference_file_name_tmp = $reference_file_name;
	$reference_file_name_tmp =~ s/.*[\\\/]//g;
	
	my $next_step = 0;
	my $repeat_v2c_check;
	
	curr_time();
	$terminal->insert('end', "Checking MD5 of the VCF file, please wait ...\n\n",);
	$terminal->see('end');
	
	$repeat_v2c_check = $mw->repeat( 100 => sub {
		if ( $v2c_check_result{err_type} ne "" )
		{
			if ( $v2c_check_result{err_type} eq "OK" )
			{
				$processing_gif->stop_animation;			
				$raw_vcf_check->configure(-image => $ok_image);
				
				curr_time();
				$terminal->insert('end', "Using previously created v2c file: '");
				$terminal->insert('end', $v2c_file_name_tmp, 'mark',);
				$terminal->insert('end', "'\n\n");
				$terminal->see('end');
				
				if ($reference_analysis_results[0] == 1 and $enzyme_analysis_results[0] == 1)
				{
					$L_lower_col1_mining_button->configure(-state => 'normal');
				}
				
				$v2c_file_name = $working_dir . $v2c_file_name_tmp;
			}
			else
			{
				if ( $v2c_check_result{err_type} eq "wrong_vcf_md5" )
				{
					curr_time();
					$terminal->insert('end', "The v2c file:  '");
					$terminal->insert('end', $v2c_file_name_tmp, 'mark',);
					$terminal->insert('end', "' in the working direcotry, is related to another VCF file. \n");
					$terminal->insert('end', "Deleting '");
					$terminal->insert('end', $v2c_file_name_tmp, 'mark',);
					$terminal->insert('end', "' ...\n\n");
					$terminal->see('end');
					
					unlink $working_dir . $v2c_file_name_tmp;
				}
				elsif ( $v2c_check_result{err_type} eq "wrong_reference_md5" )
				{
					my $reference_file_name_tmp = $reference_file_name;
					$reference_file_name_tmp =~ s/.*[\\\/]//g;
					curr_time();
					$terminal->insert('end', "The v2c file:  '");
					$terminal->insert('end', $v2c_file_name_tmp, 'mark',);
					$terminal->insert('end', "' in the working direcotry, does not match the selected reference file: '");
					$terminal->insert('end', $reference_file_name_tmp, 'mark',);
					$terminal->insert('end', "'.\n");
					$terminal->insert('end', "Deleting '");
					$terminal->insert('end', $v2c_file_name_tmp, 'mark',);
					$terminal->insert('end', "' ...\n\n");
					$terminal->see('end');
					
					unlink $working_dir . $v2c_file_name_tmp;
				}
				elsif ( $v2c_check_result{err_type} eq "no_md5" )
				{
					curr_time();
					$terminal->insert('end', "The v2c file:  '");
					$terminal->insert('end', $v2c_file_name_tmp, 'mark',);
					$terminal->insert('end', "' in the working direcotry, does not have the md5 hash header.\n");
					$terminal->insert('end', "Deleting '");
					$terminal->insert('end', $v2c_file_name_tmp, 'mark',);
					$terminal->insert('end', "' ...\n\n");
					$terminal->see('end');
					
					unlink $working_dir . $v2c_file_name_tmp;
				}
				elsif ( $v2c_check_result{err_type} eq "cannot_open_vcf" )
				{
					$processing_gif->stop_animation;			
					$raw_vcf_check->configure(-image => $warning_image);
					
					curr_time();
					$terminal->insert('end', "Cannot open the file: '");
					$terminal->insert('end', "$raw_vcf_file_name_tmp", 'mark');
					$terminal->insert('end', "'\n\n");
					$terminal->see('end');
					
					$raw_vcf_analyze_button->configure(-state => 'normal');
					$raw_chooseFile_button->configure(-state => 'normal');
					$repeat_v2c_check->cancel;
					
					warning_dial("Cannot open the file ::: $v2c_check_result{err_value}");
				}
				elsif ( $v2c_check_result{err_type} eq "cannot_open_v2c" )
				{
					$processing_gif->stop_animation;			
					$raw_vcf_check->configure(-image => $warning_image);
					
					curr_time();
					$terminal->insert('end', "Warning", 'warning');
					$terminal->insert('end', " - cannot open the file: '");
					$terminal->insert('end', "$v2c_file_name_tmp", 'mark');
					$terminal->insert('end', "'\n\n");
					$terminal->see('end');
					
					$raw_vcf_analyze_button->configure(-state => 'normal');
					$raw_chooseFile_button->configure(-state => 'normal');
					$repeat_v2c_check->cancel;
					
					warning_dial("Cannot open the file ::: $v2c_check_result{err_value}");
				}
				
				if ( $v2c_check_result{err_type} eq "cannot_open_vcf" or $v2c_check_result{err_type} eq "cannot_open_v2c" )
				{					
					return;
				}				
				
				curr_time();
				$terminal->insert('end', "Start VCF file '");
				$terminal->insert('end', "$raw_vcf_file_name_tmp", 'mark');
				$terminal->insert('end', "' integrity check ...\n\n");
				$terminal->see('end');
				
				$jobID = 3;
				my $repeat;
				$repeat = $mw->repeat( 100 => sub {
					if ($vcf_analysis_results{err_code} != 0)
					{					
						if ($vcf_analysis_results{err_code} == 1)
						{
							curr_time();
							$terminal->insert('end', "Integrity check of the VCF file '");
							$terminal->insert('end', $raw_vcf_file_name_tmp, 'mark',);
							$terminal->insert('end', "' confirmed. Statistics:\n- number of individuals: ");
							$terminal->insert('end', $vcf_analysis_results{NoOfIndv}, 'mark');
							$terminal->insert('end', "\n- number of SNPs/InDels: ",);
							$terminal->insert('end', "$vcf_analysis_results{NoOfSNPs}\n\n", 'mark');
							$terminal->see('end');
							
							curr_time();
							$terminal->insert('end', "Converting the '");
							$terminal->insert('end', $raw_vcf_file_name_tmp, 'mark');
							$terminal->insert('end', "' file into '");
							$terminal->insert('end', $v2c_file_name_tmp, 'mark');
							$terminal->insert('end', "' file ... ");
							$terminal->insert('end', "0.0%", 'percent');
							$terminal->see('end');
							$next_step = 1;
							$jobID = 4;
							@seqExtractor_error_code = ("",0);
							
							my $percent_index_start;
							
							$percent_index_start = $terminal->search(-regexp, -backwards => '[0-9]+\.[0-9]%', 'end');
							my $repeat2;
							$repeat2 = $mw->repeat( 100 => sub {
								if ($seqExtractor_error_code[1] != 0)
								{								
									my $percent = ( $line_vcf / $vcf_analysis_results{NoOfSNPs}) * 100;
									my $percent_index_stop = $terminal->search(-regexp, -backwards => '%', 'end');
									$percent_index_stop += 1;
									$terminal->delete("$percent_index_start", "$percent_index_stop");
									$terminal->insert("$percent_index_start", sprintf ("%.1f%%", $percent));
									
									$processing_gif->stop_animation;
																		
									my $sequencesNotPresentInRef_No = scalar(@sequencesNotPresentInRef);
									my $sequencesNotPresentInRef_No_forReport = scalar(@sequencesNotPresentInRef);
									my $markersOnTheEdge_No = scalar(@markersOnTheEdge);
									my $markersOnTheEdge_No_forReport = scalar(@markersOnTheEdge);
									
									if ($seqExtractor_error_code[1] == 1 and $sequencesNotPresentInRef_No == 0 and $markersOnTheEdge_No == 0 )
									{
										$raw_vcf_check->configure(-image => $ok_image);
										
										$terminal->insert('end', "\n\n");
										curr_time();
										$terminal->insert('end', "Convertion completed sucessfully.\n\n");
										$terminal->see('end');
										
										if ($reference_analysis_results[0] == 1 and $enzyme_analysis_results[0] == 1 and $vcf_analysis_results{err_code} == 1)
										{
											$L_lower_col1_mining_button->configure(-state => 'normal');
										}
										
										$v2c_file_name = $working_dir . $v2c_file_name_tmp;
									}
									elsif ($seqExtractor_error_code[1] == 1 and $snpsNo > 0 and ($sequencesNotPresentInRef_No > 0 or $markersOnTheEdge_No > 0) )
									{
										
										$raw_vcf_check->configure(-image => $ok_image);
										
										$terminal->insert('end', "\n\n");
										
										curr_time();
										$terminal->insert('end', "Warning",'warning');
										if ( $sequencesNotPresentInRef_No > 0 and $markersOnTheEdge_No > 0)
										{
											$terminal->insert('end', " - convertion completed. However, some problems occurred:\n");
										}
										else
										{
											$terminal->insert('end', " - convertion completed. However, a problem occurred:\n");
										}
										
										if ( $sequencesNotPresentInRef_No > 0 )
										{
											if ( $sequencesNotPresentInRef_No <= 5 )
											{
												$terminal->insert('end', "The SNPs/InDels listed below were located in the sequences not present in the reference file '");
												$terminal->insert('end', "$reference_file_name_tmp", 'mark');
												$terminal->insert('end', "':\n");
													
												foreach my $orphanedSeq (@sequencesNotPresentInRef)
												{
													$terminal->insert('end', "- ");
													$terminal->insert('end', "$orphanedSeq\n",'mark');
												}
												$terminal->see('end');
											}
											else
											{
												LOG("\n# Convertion of VCF file into v2c format");
												LOG("# SNPs/InDels located in the sequences not present in the reference file $reference_file_name_tmp");
												foreach my $text (@sequencesNotPresentInRef)
												{
													LOG($text);
												}
												$terminal->insert('end', "$sequencesNotPresentInRef_No_forReport", 'mark');
												$terminal->insert('end', " SNPs/InDels were located in the sequences not present in the reference file '");
												$terminal->insert('end', "$reference_file_name_tmp", 'mark');
												$terminal->insert('end', "'. The list of those SNPs/InDels were saved in the ");
												$terminal->insert('end', "log.txt", 'mark');
												$terminal->insert('end', " file.\n");
												$terminal->see('end');
											}
										}
										
										if ( $markersOnTheEdge_No > 0)
										{															
											if ( $markersOnTheEdge_No <= 5 )
											{
												$terminal->insert('end', "The SNPs/InDels listed below were located closer to the edges of the sequences than ");
												$terminal->insert('end', "$snps_seq_len", 'mark');
												$terminal->insert('end', " bp:\n");
												
												foreach my $markerOnTheEdge (@markersOnTheEdge)
												{
													$terminal->insert('end', "- ");
													$terminal->insert('end', "$markerOnTheEdge\n",'mark');
												}
												$terminal->see('end');
											}
											else
											{
												LOG("\n# Convertion of VCF file into v2c format");
												LOG("# SNPs/InDels located closer to the edges of the sequences than 40 bp");
												foreach my $text (@markersOnTheEdge)
												{
													LOG($text);
												}
												$terminal->insert('end', "$markersOnTheEdge_No_forReport", 'mark');
												$terminal->insert('end', " SNPs/InDels were located closer to the edges of the sequences than ");
												$terminal->insert('end', "$snps_seq_len", 'mark');
												$terminal->insert('end', " bp:\n");
												$terminal->insert('end', "The list of those SNPs/InDels were saved in the ");
												$terminal->insert('end', "log.txt", 'mark');
												$terminal->insert('end', " file.\n");
												$terminal->see('end');
											}
										}
										
										$terminal->insert('end', "\n");
										$terminal->see('end');
										
										if ($reference_analysis_results[0] == 1 and $enzyme_analysis_results[0] == 1 and $vcf_analysis_results{err_code} == 1)
										{
											$L_lower_col1_mining_button->configure(-state => 'normal');
										}
										
										$v2c_file_name = $working_dir . $v2c_file_name_tmp;
									}
									elsif ($seqExtractor_error_code[1] == 3)
									{
										$raw_vcf_check->configure(-image => $fail_image);
										
										$terminal->insert('end', "\n\n");
										curr_time();
										$terminal->insert('end', "Warning", 'warning');
										$terminal->insert('end', " - cannot open the index file '");
										$terminal->insert('end', "$reference_file_name_tmp.index", 'mark');
										$terminal->insert('end', "'.\n\n");
										$terminal->see('end');
										
										unlink($working_dir . $v2c_file_name_tmp) if ( -e $working_dir . $v2c_file_name_tmp );
									}
									elsif ($seqExtractor_error_code[1] == 4)
									{
										$raw_vcf_check->configure(-image => $fail_image);
										
										$terminal->insert('end', "\n\n");
										curr_time();
										$terminal->insert('end', "Warning", 'warning');
										$terminal->insert('end', " - cannot open the VCF file '");
										$terminal->insert('end', "$raw_vcf_file_name_tmp", 'mark');
										$terminal->insert('end', "'.\n\n");
										$terminal->see('end');
										
										unlink($working_dir . $v2c_file_name_tmp) if ( -e $working_dir . $v2c_file_name_tmp );
									}
									elsif ($seqExtractor_error_code[1] == 5)
									{
										$raw_vcf_check->configure(-image => $fail_image);
										
										$terminal->insert('end', "\n\n");
										curr_time();
										$terminal->insert('end', "Warning", 'warning');
										$terminal->insert('end', " - cannot write to the file '");
										$terminal->insert('end', "$v2c_file_name_tmp", 'mark');
										$terminal->insert('end', "'.\n\n");
										$terminal->see('end');
										
										unlink($working_dir . $v2c_file_name_tmp) if ( -e $working_dir . $v2c_file_name_tmp );
									}
									elsif ($snpsNo == 0)
									{
										$raw_vcf_check->configure(-image => $fail_image);
										
										$terminal->insert('end', "\n\n");
										curr_time();
										$terminal->insert('end', "Warning",'warning');
										$terminal->insert('end', " - convertion failed. Something went horribly wrong. \n");
										$terminal->insert('end', "Please, check whether you have choosen the right reference file.\n\n");
										$terminal->see('end');
										
										unlink($working_dir . $v2c_file_name_tmp) if ( -e $working_dir . $v2c_file_name_tmp );
									}
									
									$line_vcf = 0;
									$repeat2->cancel;
									
									if ( $seqExtractor_error_code[1] == 3 or $seqExtractor_error_code[1] == 4 )
									{
										warning_dial("Cannot open the file ::: $seqExtractor_error_code[0]");
									}

									if ( $seqExtractor_error_code[1] == 5 )
									{
										warning_dial("Cannot write to the file ::: $seqExtractor_error_code[0]");
									}
								}
								elsif ($line_vcf > 0 and $vcf_analysis_results{NoOfSNPs} > 0)
								{
									my $percent = ( ($line_vcf - 1) / $vcf_analysis_results{NoOfSNPs}) * 100;
									
									if (!defined $percent_index_start)
									{
										$percent_index_start = $terminal->search(-regexp, -backwards => '[0-9]+\.[0-9]%', 'end');
									}

									my $percent_index_stop = $terminal->search(-regexp, -backwards => '%', 'end');
									$percent_index_stop += 1;
									$terminal->delete("$percent_index_start", "$percent_index_stop");
									$terminal->insert("$percent_index_start", sprintf ("%.1f%%", $percent));
								}
							});
						}
						elsif ($vcf_analysis_results{err_code} == 3)
						{
							$raw_vcf_check->configure(-image => $fail_image);
							
							curr_time();
							$terminal->insert('end', "Warning", 'warning');
							$terminal->insert('end', " - the file '");
							$terminal->insert('end', "$raw_vcf_file_name_tmp", 'mark');
							$terminal->insert('end', "' does not have the header '");
							$terminal->insert('end', "#CHROM",'mark');
							$terminal->insert('end', "'. Is it really VCF file?\n\n");
							$terminal->see('end');
						}
						elsif ($vcf_analysis_results{err_code} == 6)
						{
							$raw_vcf_check->configure(-image => $warning_image);
							$L_lower_col1_mining_button->configure(-state => 'disabled');
							
							curr_time();
							$terminal->insert('end', "Warning", 'warning');
							$terminal->insert('end', " - cannot open the file: '");
							$terminal->insert('end', "$raw_vcf_file_name_tmp", 'mark');
							$terminal->insert('end', "'. Is it really VCF file?\n\n");
							$terminal->see('end');
						}
						
						$repeat->cancel;
						
						if ($vcf_analysis_results{err_code} == 6)
						{
							warning_dial("Cannot open the file ::: $vcf_analysis_results{err_value}");
						}
					}
				} );
			}
			
			$raw_vcf_analyze_button->configure(-state => 'normal');
			$raw_chooseFile_button->configure(-state => 'normal');
			$repeat_v2c_check->cancel;
		}

	} );
}


#---------------------------------------------------------------------#
# The subroutine triggers and checks the progress of caps mining step #
#---------------------------------------------------------------------#
sub start_caps_mining
{
	@caps_mining_results = ("",0);
	$jobID = 5;
	$capsMining_percent = 0;
	$enzyme_zip = 0;
	$enzyme_zip_stop = 0;
	$caps_mining_result_label->packForget;
	
	my $repeat;
	my $enzymes_for_analysis_No = keys %enzymes_for_analysis;
	
	$cfw_open_button->configure(-state => 'disabled');
	
	$repeat = $mw->repeat( 100 => sub {
		$capsMining_percent = ( ($enzyme_zip) / $enzymes_for_analysis_No) * 100;
		if ($enzyme_zip > 0 and $capsMining_percent < 100 and $caps_mining_results[1] == 0)
		{
			$caps_mining_prepare_enzymes_label->pack(-side => 'left', -anchor => 'w', -padx => 5);
			
			$caps_mining_prepare_enzymes_label->configure(-text => sprintf ("Preparing enzymes list   %.1f%%", $capsMining_percent) );
		}
		elsif ( $caps_mining_results[1] == 3 )
		{
			my $reference_file_name_tmp = $reference_file_name;
			$reference_file_name_tmp =~ s/.*[\\\/]//g;
			
			curr_time();
			$terminal->insert('end', "Warning", 'warning');
			$terminal->insert('end', " - cannot open the index file: '");
			$terminal->insert('end', "$reference_file_name_tmp", 'mark');
			$terminal->insert('end', "'.\n\n");
		
			$repeat->cancel;
			$L_lower_col1_mining_button->configure(-state => 'normal');
			$cfw_open_button->configure(-state => 'normal');
			if ( $caps_mining_results[1] == 3 )
			{
				warning_dial("Cannot open the file ::: $caps_mining_results[0]");
			}
		}
		elsif ($capsMining_percent == 100)
		{
			$caps_mining_prepare_enzymes_label->pack(-side => 'left', -anchor => 'w', -padx => 5);
			$caps_mining_prepare_enzymes_label->configure(-text => sprintf ("Preparing enzymes list   %.1f%%", $capsMining_percent) );
			$caps_mining_prepare_enzymes_label->packForget;
			$capsMining_percent = 0;
			$caps_mining_progress_label->configure(-text => sprintf ("%d/%d   %.1f%%", $actualSNPNo,$vcf_analysis_results{NoOfSNPs},$capsMining_percent) );
			$caps_mining_progress_frame->pack(-side => 'left', -anchor => 'w', -padx => 5);
			
			curr_time();
			$terminal->insert('end', "Start CAPS mining ...\n\n");
			$terminal->see('end');
			my $repeat2;
			
			$repeat2 = $mw->repeat( 100 => sub {						
				if ($caps_mining_results[1] != 0)
				{	
					@markersOnTheEdge = uniq(@markersOnTheEdge);
					
					my $markersOnTheEdge_No = scalar(@markersOnTheEdge);
					my $markersOnTheEdge_No_forReport = scalar(@markersOnTheEdge);
					
					if ($caps_mining_results[1] == 1 and $markersOnTheEdge_No == 0)
					{
						$caps_mining_progress_frame->packForget;

						$caps_mining_result_label->configure(-text => "$caps_mining_results[0] CAPS found");
						$caps_mining_result_label->pack(-side => 'left', -anchor => 'w');
						
						curr_time();						
						$terminal->insert('end', "CAPS mining finished.");
						$terminal->insert('end', " $caps_mining_results[0]", 'mark');
						$terminal->insert('end', " CAPS were saved to the file: '");
						$terminal->insert('end', "caps_markers.txt", 'mark');
						$terminal->insert('end', "'\n\n");;
						$terminal->see('end');						
					}
					elsif ($caps_mining_results[1] == 1 and $markersOnTheEdge_No > 0)
					{
						$caps_mining_progress_frame->packForget;

						$caps_mining_result_label->configure(-text => "$caps_mining_results[0] CAPS found");
						$caps_mining_result_label->pack(-side => 'left', -anchor => 'w');
												
						curr_time();
						$terminal->insert('end', "Warning",'warning');
						$terminal->insert('end', " - CAPS mining finished.");
						$terminal->insert('end', " $caps_mining_results[0]", 'mark');
						$terminal->insert('end', " CAPS were saved to the file: '");
						$terminal->insert('end', "caps_markers.txt", 'mark');
						$terminal->insert('end', "'. However, a problem occurred:\n");
												
						if ( $markersOnTheEdge_No <= 5 )
						{
							$terminal->insert('end', "The SNPs/InDels listed below were located closer to the edges of the sequences than ");
							$terminal->insert('end', "$output_seq_len", 'mark');
							$terminal->insert('end', " bp:\n");
						
							foreach my $markerOnTheEdge (@markersOnTheEdge)
							{
								$terminal->insert('end', "- ");
								$terminal->insert('end', "$markerOnTheEdge\n",'mark');
							}
							$terminal->see('end');
						}
						else
						{
							LOG("\n# CAPS mining");
							LOG("# SNPs/InDels located closer to the edges of the sequences than 500 bp");
							foreach my $text (@markersOnTheEdge)
							{
								LOG($text);
							}
							$terminal->insert('end', "$markersOnTheEdge_No_forReport", 'mark');
							$terminal->insert('end', " SNPs/InDels were located closer to the edges of the sequences than ");
							$terminal->insert('end', "$output_seq_len", 'mark');
							$terminal->insert('end', " bp:\n");
							$terminal->insert('end', "The list of those SNPs/InDels were saved to the file: ");
							$terminal->insert('end', "log.txt", 'mark');
							$terminal->insert('end', " file.\n");
							$terminal->see('end');
						}

						$terminal->insert('end', "\n");
						$terminal->see('end');						
					}
					elsif ($caps_mining_results[1] == 2)
					{
						$caps_mining_progress_frame->packForget;

						$caps_mining_result_label->configure(-text => "$caps_mining_results[0] CAPS found");
						$caps_mining_result_label->pack(-side => 'left', -anchor => 'w');
						
						curr_time();
						$terminal->insert('end', "Warning", 'warning');
						$terminal->insert('end', " - CAPS mining canceled.");
						$terminal->insert('end', " $caps_mining_results[0]", 'mark');
						$terminal->insert('end', " CAPS were saved to the file: '");
						$terminal->insert('end', "caps_markers.txt", 'mark');
						$terminal->insert('end', "'\n\n");;
						$terminal->see('end');
					}
					elsif ($caps_mining_results[1] == 2 and $markersOnTheEdge_No > 0)
					{
						$caps_mining_progress_frame->packForget;

						$caps_mining_result_label->configure(-text => "$caps_mining_results[0] CAPS found");
						$caps_mining_result_label->pack(-side => 'left', -anchor => 'w');
						
						curr_time();
						$terminal->insert('end', "Warning", 'warning');
						$terminal->insert('end', " - CAPS mining canceled.");
						$terminal->insert('end', " $caps_mining_results[0]", 'mark');
						$terminal->insert('end', " CAPS were saved to the file: '");
						$terminal->insert('end', "caps_markers.txt", 'mark');
						$terminal->insert('end', "'. However, a problem occurred:\n");
						
						if ( $markersOnTheEdge_No <= 5 )
						{
							$terminal->insert('end', "The SNPs/InDels listed below were located closer to the edges of the sequences than ");
							$terminal->insert('end', "$output_seq_len", 'mark');
							$terminal->insert('end', " bp:\n");
							
							for (my $i = 0; $i < $markersOnTheEdge_No - 1; $i++ )
							{
								$terminal->insert('end', "- ",'paragraph');
								$terminal->insert('end', "$markersOnTheEdge[$i]\n",'mark');
							}
							$terminal->insert('end', "- ");
							$terminal->insert('end', "$markersOnTheEdge[$markersOnTheEdge_No-1]\n\n",'mark');
							$terminal->see('end');
						}
						else
						{
							LOG("\n# CAPS mining");
							LOG("# SNPs/InDels located closer to the edges of the sequences than 500 bp");
							foreach my $text (@markersOnTheEdge)
							{
								LOG($text);
							}
							$terminal->insert('end', "$markersOnTheEdge_No_forReport", 'mark');
							$terminal->insert('end', " SNPs/InDels were located closer to the edges of the sequences than ");
							$terminal->insert('end', "$output_seq_len", 'mark');
							$terminal->insert('end', " bp:\n");
							$terminal->insert('end', "The list of those SNPs/InDels were saved to the file: ");
							$terminal->insert('end', "log.txt", 'mark');
							$terminal->insert('end', " file.\n");
							$terminal->see('end');
						}
												
						$terminal->insert('end', "\n");
						$terminal->see('end');
					}
					elsif ( $caps_mining_results[1] == 4 )
					{
						$caps_mining_progress_frame->packForget;
						
						my $v2c_file_name_tmp = $v2c_file_name;
						$v2c_file_name_tmp =~ s/.*[\\\/]//g;
						
						curr_time();
						$terminal->insert('end', "Warning", 'warning');
						$terminal->insert('end', " - cannot write to the output file: '");
						$terminal->insert('end', "caps_markers.txt", 'mark');
						$terminal->insert('end', "'.\n\n");
					}
					elsif ( $caps_mining_results[1] == 5 )
					{
						$caps_mining_progress_frame->packForget;
						$raw_vcf_check->configure(-image => $warning_image);
						
						my $v2c_file_name_tmp = $v2c_file_name;
						$v2c_file_name_tmp =~ s/.*[\\\/]//g;
						
						curr_time();
						$terminal->insert('end', "Warning", 'warning');
						$terminal->insert('end', " - cannot open the v2c file: '");
						$terminal->insert('end', "$v2c_file_name_tmp", 'mark');
						$terminal->insert('end', "'.\n\n");
					}
										
					$L_lower_col1_mining_button->configure(-state => 'normal');
					$cfw_open_button->configure(-state => 'normal');
					
					$repeat2->cancel;
					
					if ( $caps_mining_results[1] == 4 )
					{
						warning_dial("Cannot write to the file ::: $caps_mining_results[0]");
					}
					elsif ( $caps_mining_results[1] == 5 )
					{
						warning_dial("Cannot open the file ::: $caps_mining_results[0]");
					}
				}
				elsif ($actualSNPNo > 0 and $vcf_analysis_results{NoOfSNPs} > 0)
				{
					$capsMining_percent = ( ($actualSNPNo) / $vcf_analysis_results{NoOfSNPs}) * 100;
					
					$caps_mining_progress_label->configure(-text => sprintf ("%d/%d   %.1f%%", $actualSNPNo,$vcf_analysis_results{NoOfSNPs},$capsMining_percent) );					

				}
			});
			
			$repeat->cancel;
		}		
	});
}


#-----------------------------------------------------------------------------------#
# The subroutine triggers and checks the progress of the single-cut filtration step #
#-----------------------------------------------------------------------------------#
sub start_singleCut_filter
{
	
	$cfw_scf_CAPS_filtering_percent = 0;
	$numberOfSNPsBefore = 0;
	$cfw_scf_error_label->packForget;
	$cfw_scf_progress_label->configure(-text => sprintf ("%d/%d   %.1f%%", $numberOfSNPsBefore,$total_caps_number,$cfw_scf_CAPS_filtering_percent) );
	$cfw_scf_progress_frame->pack(-side => 'left', -anchor => 'w', -padx => 5);
	my $repeat;
	@singleCutSite_results = ("", 0);
	$jobID = 6;
	
	curr_time();
	$terminal->insert('end', "Start single-cut CAPS filtering ...\n\n");
	$terminal->see('end');
	$repeat = $mw->repeat( 100 => sub {
		if ($singleCutSite_results[1] != 0)
		{
			my $cfw_scf_output_file_tmp = $cfw_scf_output_file;
			$cfw_scf_output_file_tmp =~ s/.+\/(.+\.txt)/$1/;
				
			if ($singleCutSite_results[1] == 1)
			{
				$cfw_scf_CAPS_filtering_percent = ( $numberOfSNPsBefore / $total_caps_number ) * 100;
				$cfw_scf_progress_label->configure(-text => sprintf ("%d/%d   %.1f%%", $numberOfSNPsBefore,$total_caps_number,$cfw_scf_CAPS_filtering_percent) );

				$cfw_scf_progress_frame->packForget;
				$numberOfSNPsAfter =~ s/(\d{1,3}?)(?=(\d{3})+$)/$1 /g;
				$cfw_scf_error_label->configure(-text => "$numberOfSNPsAfter CAPS with single-cut site were saved to the file: '$cfw_scf_output_file_tmp'");
				$cfw_scf_error_label->pack(-side => 'left', -anchor => 'w');
				
				curr_time();
				$terminal->insert('end', "Single-cut CAPS filtering finished.");
				$terminal->insert('end', " $numberOfSNPsAfter", 'mark');
				$terminal->insert('end', " CAPS were saved to the file: '");
				$terminal->insert('end', "$cfw_scf_output_file_tmp", 'mark');
				$terminal->insert('end', "' file\n\n");
				$terminal->see('end');
			}
			elsif ($singleCutSite_results[1] == 2)
			{
				$cfw_scf_progress_frame->packForget;
				$numberOfSNPsAfter =~ s/(\d{1,3}?)(?=(\d{3})+$)/$1 /g;
				$cfw_scf_error_label->configure(-text => "$numberOfSNPsAfter CAPS with single-cut site were saved to the file: '$cfw_scf_output_file_tmp'");
				$cfw_scf_error_label->pack(-side => 'left', -anchor => 'w');
				
				curr_time();
				$terminal->insert('end', "Warning",'warning');
				$terminal->insert('end', " - single-cut CAPS filtering canceled.");
				$terminal->insert('end', " $numberOfSNPsAfter", 'mark');
				$terminal->insert('end', " CAPS were saved to the file: '");
				$terminal->insert('end', "$cfw_scf_output_file_tmp", 'mark');
				$terminal->insert('end', "' file\n\n");
				$terminal->see('end');
			}
			elsif ($singleCutSite_results[1] == 3)
			{
				$cfw_scf_progress_frame->packForget;
				
				curr_time();
				$terminal->insert('end', "Warning",'warning');
				$terminal->insert('end', " - cannot write to the output file: ");
				$terminal->insert('end', "$cfw_scf_output_file_tmp", 'mark');
				$terminal->insert('end', "'.\n\n");
				$terminal->see('end');
			}
			elsif ($singleCutSite_results[1] == 4)
			{
				$cfw_scf_progress_frame->packForget;
				
				my $cfw_scf_input_file_tmp = $cfw_scf_input_file;
				$cfw_scf_input_file_tmp =~ s/.*[\\\/]//g;
				$cfw_scf_inputFile_check->configure(-image => $warning_image);
				
				curr_time();
				$terminal->insert('end', "Warning",'warning');
				$terminal->insert('end', " - cannot open the file: ");
				$terminal->insert('end', "$cfw_scf_input_file_tmp", 'mark');
				$terminal->insert('end', "'.\n\n");
				$terminal->see('end');
			}
			$cfw_scf_start_button->configure(-state => 'normal');			
			
			$repeat->cancel;
			
			if ($singleCutSite_results[1] == 3)
			{
				warning_dial("Cannot write to the file ::: $singleCutSite_results[0]");
			}
			elsif ($singleCutSite_results[1] == 4)
			{
				warning_dial("Cannot open the file ::: $singleCutSite_results[0]");
				$Caps_Filtration_Window->focusForce;
			}
		}
		
		elsif ($numberOfSNPsBefore > 0 and $total_caps_number > 0)
		{
			$cfw_scf_CAPS_filtering_percent = ( $numberOfSNPsBefore / $total_caps_number ) * 100;
			$cfw_scf_progress_label->configure(-text => sprintf ("%d/%d   %.1f%%", $numberOfSNPsBefore,$total_caps_number,$cfw_scf_CAPS_filtering_percent) );			
		}
	});
}


#------------------------------------------------------------#
# The subroutine retrieves a path to the specific input file #
#------------------------------------------------------------#
sub fileDialog {
	my $w = shift;
	my $file_type = shift;
	my $types;
	my $file;
	my @types_reference =
		(["Fasta files", [qw/.fa .fasta .fna/],'TEXT'],
		["Text files", '.txt','TEXT'],
		["All files", '*']
	);
	my @types_vcf =
		(["VCF files", '.vcf', 'TEXT'],
		["Text files", '.txt'],
		["All files", '*']
	);

	if ($file_type eq "enzyme")
	{
		$file = $w->getOpenFile(-defaultextension => '.txt');
		if (defined $file and $file ne "")
		{
			$enzyme_check->configure(-image => '', -text => 'none');
			
			$enzyme_entry->delete(0, 'end');
			$enzyme_entry->insert(0, $file);
			$enzyme_entry->xview('end');
		}
	}
	elsif ($file_type eq "reference")
	{
		$file = $w->getOpenFile(-filetypes => \@types_reference, -defaultextension => '.fa');
		if (defined $file and $file ne "")
		{
			$reference_check->configure(-image => '', -text => 'none');
			$raw_vcf_check->configure(-image => '', -text => 'none');
			$vcf_analysis_results{err_code} = 0;
			$L_lower_col1_mining_button->configure(-state => 'disabled');
			
			@reference_analysis_results = (0);
			
			$reference_entry->delete(0, 'end');
			$reference_entry->insert(0, $file);
			$reference_entry->xview('end');
		}
	}
	elsif ($file_type eq "raw_vcf")
	{
		$file = $w->getOpenFile(-filetypes => \@types_vcf, -defaultextension => '.vcf');
		if (defined $file and $file ne "")
		{
			$raw_vcf_check->configure(-image => '', -text => 'none');
			$L_lower_col1_mining_button->configure(-state => 'disabled');
			
			$raw_vcf_entry->delete(0, 'end');
			$raw_vcf_entry->insert(0, $file);
			$raw_vcf_entry->xview('end');
		}
	}
	elsif ($file_type eq "vcf2caps_output")
	{
		$file = $w->getOpenFile(-defaultextension => '.txt');
		if (defined $file and $file ne "")
		{
			$cfw_gf_inputFile_check->configure(-image => '', -text => 'none', -foreground => 'grey');
			
			$cfw_gf_inputFile_entry->delete(0, 'end');
			$cfw_gf_inputFile_entry->insert(0, $file);
			$cfw_gf_inputFile_entry->xview('end');
			$cfw_gf_start_button->configure(-state => 'disabled');
			$Caps_Filtration_Window->focusForce;
		}
	}
	elsif ($file_type eq "scf_vcf2caps_output")
	{
		$file = $w->getOpenFile(-defaultextension => '.txt');
		if (defined $file and $file ne "")
		{
			$cfw_scf_inputFile_check->configure(-image => '', -text => 'none', -foreground => 'grey');
			
			$cfw_scf_inputFile_entry->delete(0, 'end');
			$cfw_scf_inputFile_entry->insert(0, $file);
			$cfw_scf_inputFile_entry->xview('end');
			$cfw_scf_start_button->configure(-state => 'disabled');
			$Caps_Filtration_Window->focusForce;
		}
	}
	elsif ($file_type eq "c2f_vcf2caps_output")
	{
		$file = $w->getOpenFile(-defaultextension => '.txt');
		if (defined $file and $file ne "")
		{
			$cfw_c2f_inputFile_check->configure(-image => '', -text => 'none', -foreground => 'grey');
			
			$cfw_c2f_inputFile_entry->delete(0, 'end');
			$cfw_c2f_inputFile_entry->insert(0, $file);
			$cfw_c2f_inputFile_entry->xview('end');
			$cfw_c2f_start_button->configure(-state => 'disabled');
			$Caps_Filtration_Window->focusForce;
		}
	}
}

#------------------------------------------------------------------------------------------------------------------------#
# The main function that works in the background performing specific tasks depending on the value of the variable $jobID #
#------------------------------------------------------------------------------------------------------------------------#
sub work
{
	while(1)
	{
		if ($jobID == 1)
		{
			if (defined $reference_file_name)
			{
				my $genome_exists = 0;
				my @genome_error = (0);
				if (-e $reference_file_name)
				{
					$genome_exists = 1;
					my $reference_file_name_tmp = $reference_file_name;
					$reference_file_name_tmp =~ s/\..+$//;
					$reference_file_name_tmp =~ s/.*[\\\/]//g;
					my $line_len = 0;
					my $chrom_ID = "";
					my $chrom_len = 0;
					my $last_line_len = 0;
					my $alert = 0;
					my @all_chrom_ID = ();
					
					my $fh;
					if ( !open $fh, "<", $reference_file_name )
					{
						$reference_analysis_results[1] = $!;
						$reference_analysis_results[0] = 5;
						$jobID = 0;
						next;
					}
					
					my $Ofh;
					if ( !open $Ofh, '>>', $working_dir . $reference_file_name_tmp . ".index" )
					{
						$reference_analysis_results[1] = $!;
						$reference_analysis_results[0] = 6;
						$jobID = 0;
						next;
					}
					
						print $Ofh "#$reference_md5\n";
						
						L: while (<$fh>)
						{
							chomp $_;
							if ($_ =~ /^>/)
							{
								if ($chrom_len != 0) {print $Ofh "$last_line_len\t$chrom_len\n"; $chrom_len = 0}
								$chrom_ID = $_;
								if (scalar(@all_chrom_ID) > 0)
								{
									foreach my $c (@all_chrom_ID)
									{
										if ($chrom_ID eq $c)
										{
											close $Ofh;
											unlink $working_dir . $reference_file_name_tmp . ".index";

											@genome_error = (1,$chrom_ID);
											last L;
										}
									}
								}
								push @all_chrom_ID, $chrom_ID;
								print $Ofh $_ . "\t" . tell($fh) . "\t";
								$line_len = 0;
								$alert = 0;
							}
							else
							{
								if ($line_len == 0)
								{
									$line_len = (split("", $_));
									$last_line_len = $line_len;
									$chrom_len += $line_len;
									print $Ofh "$line_len\t";
								}
								else
								{
									if ((split("", $_)) > 0)
									{
										if ($line_len != (split("", $_)))
										{
											if ($alert == 1)
											{
												close $Ofh;
												unlink $working_dir . $reference_file_name_tmp . ".index";
												
												@genome_error = (2,"$chrom_ID,$.");
												last L;
											}
											$alert = 1;
										}
										else
										{
											if ($alert == 1)
											{
												close $Ofh;
												unlink $working_dir . $reference_file_name_tmp . ".index";
												
												@genome_error = (2,"$chrom_ID,$.");
												last L;
											}
										}
										$chrom_len += (split("", $_));
										$last_line_len = (split("", $_));
									}
								}
							}
						}
						if ($chrom_len != 0) {print $Ofh "$last_line_len\t$chrom_len\n"; $chrom_len = 0}

					close $fh;
					close $Ofh;
					
				} else {$genome_exists = 0}

				if ($genome_exists == 1 and $genome_error[0] == 0) { $reference_analysis_results[0] = 1 }
				elsif ($genome_exists == 0) { $reference_analysis_results[0] = 2 }
				elsif ($genome_exists == 1 and $genome_error[0] == 1) { @reference_analysis_results = (3,$genome_error[1]) }
				elsif ($genome_exists == 1 and $genome_error[0] == 2) { @reference_analysis_results = (4,$genome_error[1]) }
			}

			$jobID = 0;
		}
		elsif ($jobID == 2)
		{
			my $date_check = 0;
			my @enzymes_tmp;
			my @enzymes_uniq_seq;
			@allEnzymesNames = ();

			if (defined $enzyme_file_name)
			{
				my $fh;
				my $line = 0;
				my $enzymes_exists = 1;
				my $db_OK = 0;
				$enzymes_db{date} = "unknown";
				my $enzymes_begin = 0;
				
				if ( !open $fh, '<', $enzyme_file_name )
				{
					$enzyme_analysis_results[1] = $!;
					$enzyme_analysis_results[0] = 4;
					$jobID = 0;
					next;
				}
								
					my @id_vs_companyName = (); # <company_ID>,<company_name>
										
					while (<$fh>)
					{
						if ($_ =~ /REBASE\sversion/ and $db_OK == 0)
						{
							$db_OK = 1;
						}
						elsif ($db_OK == 0)
						{
							next;
						}
						
						chomp $_;						
						my @data = split(" ", $_);
						
						if (defined $data[0])
						{
							if ($date_check == 0)
							{
								for (my $i = 0; $i < scalar(@data); $i++)
								{
									foreach my $month (qw/Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec/)
									{
										if ($month eq $data[$i])
										{
											$enzymes_db{date} = "$data[$i] $data[$i + 1] $data[$i + 2]";
											$date_check = 1;
											last;
										}
									}
								}
							}
							
							if ($data[0] =~ /^\.{2}$/)
							{
								$enzymes_begin = 1;
								$date_check = 1;
							}
							elsif ($data[0] =~ /^[A-Z]$/)
							{
								my $ID = $data[0];
								my $company_name = "";
								for ( my $i = 1; $i < scalar(@data); $i++ )
								{
									$company_name = $company_name . " " . $data[$i];
								}
								$company_name =~ s/ \([0-9]{1,2}\/[0-9]{1,2}\)//;
								$company_name =~ s/^ //;
								my $tmp = join(",", $ID,$company_name);
								push @id_vs_companyName, $tmp;
								$enzymes_db{companies} = join("\t",@id_vs_companyName);
							}
							elsif ($data[0] =~ /[A-Z][a-z]{2}/ and $enzymes_begin == 1)
							{
								$data[0] =~ s/;//;
								my $raw_seq = $data[2];
								$data[2] =~ s/[_']//g;
								my @tmp = (); # <cut_site>,<sequence>,<overhang>,<isoschisomers>,<company_ID>
								$data[2] =~ s/(^[nN]+)|([nN]+$)//g; # removing non-specific nucleotides (N) from beginning and end of recognized sequence
								@tmp = ($data[1],$data[2],$data[3]);
								if ($data[5] =~ /^>/)
								{
									$data[5] =~ s/>//;
									if ($data[5] eq "") { $data[5] = "null" }
									push @tmp, ("null", $data[5]);
								}
								elsif ($data[5] eq "!")
								{
									if ($data[6] =~ /^>/)
									{
										$data[6] =~ s/>//;
										if ($data[6] eq "") { $data[6] = "null" }
										push @tmp, ("null", $data[6]);
									}
									else
									{
										$data[7] =~ s/>//;
										if ($data[7] eq "") { $data[7] = "null" }
										push @tmp, ($data[6],$data[7]);
									}
								}
								else
								{
									$data[6] =~ s/>//;
									if ($data[6] eq "") { $data[6] = "null" }
									push @tmp, ($data[5],$data[6]);
								}
								
								$enzymes_db{$data[0]} = join("\t",@tmp,$raw_seq);
								
								$data[2] = uc $data[2];
								push @enzymes_tmp, join("\t", $data[2],$data[0]);
							}							
						}
					}
				close $fh;
				
				if ($db_OK == 1)
				{
					foreach my $enzyme (@enzymes_tmp)
					{
						push @enzymes_uniq_seq, ( split("\t", $enzyme) )[0];
					}
					@enzymes_uniq_seq = uniq(@enzymes_uniq_seq);
					
					foreach my $enzyme (@enzymes_tmp)
					{
						push @allEnzymesNames, split(",", ( split("\t", $enzyme) )[1] );
					}
					@allEnzymesNames = uniq(@allEnzymesNames);
					
					foreach my $enzyme_seq (@enzymes_uniq_seq)
					{
						my @data = ();
						my $seq_len = 0;
						for ( my $y = 0; $y < scalar(@enzymes_tmp); $y++ )
						{
							if ( $enzyme_seq eq ( split("\t",$enzymes_tmp[$y]) )[0] )
							{
								push @data, ( split("\t",$enzymes_tmp[$y]) )[1];
								$seq_len = scalar( split("", ( split("\t",$enzymes_tmp[$y]) )[0] ) );
							}
						}
						
						$enzymes_for_analysis{join(",",@data)} = join("\t",$enzyme_seq,$seq_len);

					}			
				}
				if ($enzymes_exists == 1 and $db_OK == 1) { @allEnzymesNames = sort{$a cmp $b} @allEnzymesNames; $enzyme_analysis_results[0] = 1 }
				elsif ($enzymes_exists == 1 and $db_OK == 0) { $enzyme_analysis_results[0] = 3 }
			}
			
			$jobID = 0;
		}
		elsif ($jobID == 3)
		{
			my $vcf_number_of_variants = 0;
			my $vcf_error = 0;
			my $vcf_exists = 0;
			@linie = ();
			my $number_of_individuals = 0;
			
			if (-e $raw_vcf_file_name)
			{
				$vcf_exists = 1;

				my $check_vcf;
				if ( !open $check_vcf, '<', "$raw_vcf_file_name" )
				{
					%vcf_analysis_results = ( err_code => 6, err_value => $! );
					
					$jobID = 0;
					next;
				}
				
				my $begin = 0;
				
				while (my $check_bier = <$check_vcf>)
				{
					chomp $check_bier;
					my @data = split("\t", $check_bier);

					if ($data[0] =~ /^#CHROM/)
					{
						$begin = 1;
						
						for (my $i = 9; $i < scalar(@data); $i++)
						{
							push @linie, $data[$i];
						}
						
						$number_of_individuals = @linie;
					}
					elsif ( $begin == 1 )
					{
						$vcf_number_of_variants++;
					}
				}
				
				if ($begin == 0)
				{
					$vcf_error = 1;
				}					
					
				close $check_vcf;
			} else { $vcf_exists = 0 }
			
			if ($vcf_exists == 1 and $vcf_error == 0) { %vcf_analysis_results = (err_code => 1, NoOfIndv => $number_of_individuals, NoOfSNPs => $vcf_number_of_variants) }
			elsif ($vcf_exists == 0) { $vcf_analysis_results{err_code} = 2 }
			elsif ($vcf_exists == 1 and $vcf_error == 1) { $vcf_analysis_results{err_code} = 3 }

			$jobID = 0;
		}
		elsif ($jobID == 4)
		{
			@seqExtractor_error_code = (vcf2snps($snps_seq_len,$reference_file_name));
			
			$jobID = 0;
		}
		elsif ($jobID == 5)
		{
			@caps_mining_results = ( caps_miner($output_seq_len,$reference_file_name) );
			$jobID = 0;
		}
		elsif ($jobID == 6)
		{
			@singleCutSite_results = ( singleCutSite() );
			
			$jobID = 0;
		}
		elsif ($jobID == 7)
		{
			$vcf_md5 = "";
			@linie = ();
			my $v2c_file_name_tmp = $raw_vcf_file_name;
			$v2c_file_name_tmp =~ s/\..+$/\.v2c/;
			$v2c_file_name_tmp =~ s/.*[\\\/]//g;
			
			my $number_of_individuals = 0;
			my $number_of_variants = 0;
			
			if ( -e $raw_vcf_file_name )
			{
				my $raw_vcf_file;
				if ( !open $raw_vcf_file, '<', $raw_vcf_file_name )
				{
					$v2c_check_result{err_value} = $!;
					$v2c_check_result{err_type} = "cannot_open_vcf";
					
					$jobID = 0;
					next;
				}
				
				binmode($raw_vcf_file);
				my $md5 = Digest::MD5->new;
				$md5->addfile($raw_vcf_file);
				$vcf_md5 = $md5->hexdigest;
				
				close $raw_vcf_file;
			}

			my $md5_ok = 0;
			if ( -e $working_dir . $v2c_file_name_tmp )
			{
				my $fh;
				if ( !open $fh, '<', $working_dir . $v2c_file_name_tmp )
				{
					$v2c_check_result{err_value} = $!;
					$v2c_check_result{err_type} = "cannot_open_v2c";
				
					$jobID = 0;
					next;
				}

				my $line_number = 1;
				while (my $line = <$fh>)
				{
					chomp $line;
					if ( $line =~ /^#/ )
					{
						$line =~ s/^#//;
						
						if ( $line_number == 1 )
						{
							if ( $line eq $reference_md5 )
							{
								$md5_ok = 1;			
							}
							else
							{
								$v2c_check_result{err_type} = "wrong_reference_md5";						
								last;
							}
						}
						elsif ( $line_number == 2 )
						{													
							if ( $line eq $vcf_md5 )
							{
								$md5_ok = 1;
							}
							else
							{
								$v2c_check_result{err_type} = "wrong_vcf_md5";								
								last;
							}
						}						
					}
					elsif ( $line =~ /^>/ )
					{
						$number_of_variants++;						
					}
					
					$line_number++;
				}
				close $fh;
			}
			else
			{
				$v2c_check_result{err_type} = "not_exist";
			}
			
			if ( $md5_ok == 1 )
			{				
				my $fh;
				if ( !open $fh, '<', $raw_vcf_file_name )
				{
					$v2c_check_result{err_value} = $!;
					$v2c_check_result{err_type} = "cannot_open_vcf";

					$jobID = 0;
					next;
				}
				
				while ( my $line = <$fh> )
				{
					chomp $line;
					if ( $line =~ /^#CHROM/ )
					{
						my @data = split("\t", $line);
						
						for (my $i = 9; $i < scalar(@data); $i++)
						{
							push @linie, $data[$i];
						}
						
						$number_of_individuals = @linie;						
						last;
					}
				}
				close $fh;
				
				$v2c_check_result{err_type} = "OK";
				
				$vcf_analysis_results{err_code} = 1;
				$vcf_analysis_results{NoOfIndv} = $number_of_individuals;
				$vcf_analysis_results{NoOfSNPs} = $number_of_variants;
			}
			elsif ( $v2c_check_result{err_type} eq "" )
			{
				$v2c_check_result{err_type} = "no_md5";
			}
						
			$jobID = 0;
		}
		elsif ($jobID == 8)
		{
			$cfw_gf_output_file = $cfw_gf_input_file;
			$cfw_gf_output_file =~ s/\.txt$/_gf\.txt/;
			$cfw_gf_output_file =~ s/.*[\\\/]//g;
			$cfw_gf_output_file = $working_dir . $cfw_gf_output_file;
			
			if (-e $cfw_gf_output_file) { unlink $cfw_gf_output_file }
					
			$stop = 0;
			
			my %cfw_groups_thread = ();
			my %data_filtered = ();
			
			for ( my $i = 1; $i < 4; $i++ )
			{
				if ( $cfw_groups{$i}{indv} )
				{
					$cfw_groups_thread{$i}{max_err} = $cfw_groups{$i}{max_err} / 100;
					$cfw_groups_thread{$i}{indv} = [ split(/[, \t\n]+/, $cfw_groups{$i}{indv}) ];
				}			
			}
			
			my $Ofh;
			if ( !open $Ofh, '>>', $cfw_gf_output_file )
			{
				@caps_filtration_result = ( 2, $! );
				
				$jobID = 0;
				next;
			}
			
			print $Ofh "#VCF2CAPS_output_file\n\n";
			close $Ofh;
			
			my $fh;
			if ( !open $fh, '<', $cfw_gf_input_file )
			{
				@caps_filtration_result = ( 3, $! );
				
				$jobID = 0;
				next;
			}
			
			my @write2file_filtered_CAPS_result = ();
			while (<$fh>)
			{
				if ( $stop == 0 )
				{
					chomp $_;
					my @data = split("\t", $_);
					if (! defined $data[0]) {next}
					if ($data[0] =~ /^>/)
					{
						$caps_filtered++;
						if (keys %data_filtered)
						{
							my $code = check_genotypes_vs_filters(\%cfw_groups_thread, \%data_filtered);
							if ($code == 1)
							{
								@write2file_filtered_CAPS_result = write2file_filtered_CAPS(\%data_filtered);
								if ( $write2file_filtered_CAPS_result[1] == 6 )
								{
									@caps_filtration_result = ( 2, $write2file_filtered_CAPS_result[0] );
									
									$jobID = 0;
									next;
								}
								
								$caps_filtration_result[1]++;
							}
						}

						%data_filtered = ();
						$data_filtered{id} = $data[0];
					}
					elsif ($data[0] =~ /^[A-Z][a-z]{2}/)
					{
						$data_filtered{enz} = $_;
					}
					elsif ($data[0] =~ /^ref/)
					{
						$data_filtered{ref} = $_;
					}
					elsif ($data[0] =~ /^alt/)
					{
						if ($data_filtered{alt})
						{
							$data_filtered{alt} = $data_filtered{alt} . ":::" . $_;
						}
						else
						{
							$data_filtered{alt} = $_;
						}
					}
					elsif ($data[0] =~ /^[0-9\.]/)
					{
						my @spliced = @data;
						splice @spliced, 0, 2;
						
						if ($data_filtered{genot_print})
						{
							$data_filtered{genot_print} = $data_filtered{genot_print} . ":::" . $_;
						}
						else
						{
							$data_filtered{genot_print} = $_;
						}
						
						if ($data[0] =~ /^[0-9]/)
						{
							if (exists $data_filtered{genot}{$data[1]})
							{
								$data_filtered{genot}{$data[1]} = $data_filtered{genot}{$data[1]} . "\t" . join("\t", @spliced);
							}
							else
							{
								$data_filtered{genot}{$data[1]} = join("\t", @spliced);
							}
						}
					}
				}
			}
			
			if (keys %data_filtered)
			{
				my $code = check_genotypes_vs_filters(\%cfw_groups_thread, \%data_filtered);
				if ($code == 1)
				{
					write2file_filtered_CAPS(\%data_filtered);
					$caps_filtration_result[1]++;
				}
			}
			close $fh;
			
			$caps_filtration_result[0] = 1;
			$jobID = 0;		
		}
		elsif ($jobID == 9)
		{
			$total_caps_number = 0;
			
			my $fh;
			
			if ( !open $fh, '<', $cfw_gf_input_file )
			{
				$vcf2capsOutput_results{err_value} = $!;
				$vcf2capsOutput_results{err_code} = 6;
				
				$jobID = 0;
				next;
			}
			
			while (<$fh>)
			{
				chomp $_;
				if ( $_ eq "#VCF2CAPS_output_file" )
				{
					$vcf2capsOutput_results{err_code} = 1;
				}
				elsif ( $_ =~ /^>/ )
				{
					$total_caps_number++;
				}
			}
			if ( $vcf2capsOutput_results{err_code} == 0 )
			{
				$vcf2capsOutput_results{err_code} = 2;
			}
			
			close $fh;
			$jobID = 0;
		}
		elsif ($jobID == 10)
		{
			$total_caps_number = 0;
			my $file_ok = 0;
			
			my $fh;
			
			if ( !open $fh, '<', $cfw_scf_input_file )
			{
				$vcf2capsOutput_results{err_value} = $!;
				$vcf2capsOutput_results{err_code} = 6;
				
				$jobID = 0;
				next;
			}
			
			while (<$fh>)
			{
				chomp $_;
				if ( $_ eq "#VCF2CAPS_output_file" )
				{
					$file_ok = 1;
				}
				elsif ( $_ =~ /^>/ )
				{
					$total_caps_number++;
				}
			}
			
			if ( $file_ok == 1 )
			{
				$vcf2capsOutput_results{err_code} = 1;
			}
			
			if ( $vcf2capsOutput_results{err_code} == 0 )
			{
				$vcf2capsOutput_results{err_code} = 2;
			}
			
			close $fh;
			$jobID = 0;
		}
		elsif ($jobID == 11)
		{
			$total_caps_number = 0;
			my %variants_name_uniq = ();
			my $file_ok = 0;
			
			my $fh;
			if ( !open $fh, '<', $cfw_c2f_input_file )
			{
				$vcf2capsOutput_results{err_value} = $!;
				$vcf2capsOutput_results{err_code} = 6;
				
				$jobID = 0;
				next;
			}
			
			while (<$fh>)
			{
				chomp $_;
				if ( $_ eq "#VCF2CAPS_output_file" )
				{
					$file_ok = 1;
				}
				elsif ( $_ =~ /^>/ and !exists $variants_name_uniq{$_} )
				{
					$total_caps_number++;
					$variants_name_uniq{$_} = 1;
				}
			}
			
			if ( $file_ok == 1 )
			{
				$vcf2capsOutput_results{err_code} = 1;
			}
			
			if ( $vcf2capsOutput_results{err_code} == 0 )
			{
				$vcf2capsOutput_results{err_code} = 2;
			}
			
			close $fh;
			$jobID = 0;
		}
		elsif ($jobID == 12)
		{
			$cfw_c2f_output_file = $cfw_c2f_input_file;
			$cfw_c2f_output_file =~ s/\.txt$/\.fasta/;
			$cfw_c2f_output_file =~ s/.*[\\\/]//g;
			$cfw_c2f_output_file = $working_dir . $cfw_c2f_output_file;
			
			if (-e $cfw_c2f_output_file) { unlink $cfw_c2f_output_file }
			
			$stop = 0;
			my %variants_name_uniq = ();
			my $variant = "";
			
			my $Ofh;
			if ( !open $Ofh, '>>', $cfw_c2f_output_file )
			{
				$caps_to_fasta_result[0] = 2;
				$caps_to_fasta_result[1] = $!;
				
				$jobID = 0;
				next;
			}
			
			my $fh;
			if ( !open $fh, '<', $cfw_c2f_input_file )
			{
				$caps_to_fasta_result[0] = 3;
				$caps_to_fasta_result[1] = $!;
				
				$jobID = 0;
				next;
			}
				
			while (<$fh>)
			{
				if ($stop == 0)
				{
					chomp $_;
				
					if ( $_ =~ /^>/ and !exists $variants_name_uniq{$_} )
					{
						print $Ofh "$_\n";
						$variants_name_uniq{$_} = 1;
						$variant = $_;
						$caps_filtered++;
						$caps_to_fasta_result[1]++;
					}
					elsif ( $_ =~ /^ref/ and $variant ne "" )
					{
						print $Ofh ( split("\t", $_) )[4] . "\n";
						$variant = "";
					}
				}
			}
			
			close $fh;
			close $Ofh;
			$caps_to_fasta_result[0] = 1;
			$jobID = 0;
		}
		elsif ($jobID == 13)
		{
			$reference_md5 = "";
			
			if ( -e $reference_file_name )
			{
				my $reference_file;
				
				if ( !open $reference_file, '<', $reference_file_name )
				{
					$reference_index_check_result{error_value} = $!;
					$reference_index_check_result{error_code} = "cannot_open_ref";
					
					$jobID = 0;
					next;
				}
				
				binmode($reference_file);
				my $md5 = Digest::MD5->new;
				$md5->addfile($reference_file);
				$reference_md5 = $md5->hexdigest;
				close $reference_file;
			}
			
			my $reference_index_file_name = $reference_file_name;
			$reference_index_file_name =~ s/\..+$/\.index/g;
			$reference_index_file_name =~ s/.*[\\\/]//g;
			$reference_index_file_name = $working_dir . $reference_index_file_name;
			
			if ( -e $reference_index_file_name )
			{
				my $fh;

				if ( !open $fh, '<', $reference_index_file_name )
				{
					$reference_index_check_result{error_value} = $!;
					$reference_index_check_result{error_code} = "cannot_open_index";
					
					$jobID = 0;
					next;
				}
				
				while (my $line = <$fh>)
				{
					chomp $line;
					
					if ( $line =~ /^#/ )
					{
						$line =~ s/^#//;
						
						if ( $line eq $reference_md5 )
						{
							$reference_index_check_result{error_code} = "OK";					
						}
						else
						{
							$reference_index_check_result{error_code} = "wrong_md5";
						}
						
						last;
					}
				}
				close $fh;
			}
			else
			{
				$reference_index_check_result{error_code} = "not_exist";
			}
			
			if ( $reference_index_check_result{error_code} eq "OK" )
			{
				$reference_analysis_results[0] = 1;
			}
			elsif ( $reference_index_check_result{error_code} eq "" )
			{
				$reference_index_check_result{error_code} = "no_md5";
			}
			
			$jobID = 0;
		}
		elsif ($jobID == 14)
		{
			my $url = 'http://rebase.neb.com/rebase/link_gcg';
			my $file = $working_dir . "link_gcg";

			my $response_code = getstore($url, $file);
			if ($response_code != 200)
			{
				$download_enzyme_db_result = 2;
			}
			elsif ($response_code == 200)
			{
				$enzyme_file_name = $file;				
				$download_enzyme_db_result = 1;
			}
			
			$jobID = 0;
		}
		sleep 1;
	}
}


#-------------------------------------------------#
# The subroutine performs filtration by genotype  #
#-------------------------------------------------#
sub check_genotypes_vs_filters
{
	(my $cfw_groups_thread, my $data_filtered) = @_;
	my $num_filters = scalar(keys %$cfw_groups_thread);
	my %filters_results = ();
	my %filters_matched_indv = ();
	my $err_No = 0;
	my $number_of_matched_genotypes = 0;
	
	my $number_of_genotypes = scalar( keys %{$$data_filtered{genot}} );
	my %used_filteres = ();
	
	GENOT: foreach my $genot ( grep { $_ =~ /^\[/ } keys %{$$data_filtered{genot}} )
	{
		my $duplicate = 0;
		my %indviduals_per_genot = map { $_ => 1 } split("\t", $$data_filtered{genot}{$genot});
		
		foreach my $filter ( sort keys %$cfw_groups_thread )
		{
			$err_No = 0;

			my $No_of_indv_filter_genotypes = scalar( @{$$cfw_groups_thread{$filter}{indv}} );
			
			foreach my $indv ( @{$$cfw_groups_thread{$filter}{indv}} )
			{
				if ( ! exists($indviduals_per_genot{$indv}) )
				{
					$err_No++;
				}
				else
				{
					if ( !defined $filters_matched_indv{$genot}{$filter} )
					{
						$filters_matched_indv{$genot}{$filter} = 1;
					}
					else
					{
						$filters_matched_indv{$genot}{$filter}++;
					}
					
				}
			}
			
			my $err_percent = $err_No / $No_of_indv_filter_genotypes;

			if ( $err_percent <= $$cfw_groups_thread{$filter}{max_err} and !exists( $used_filteres{$filter}  ) )
			{
				$filters_results{$genot} = $filter;
				$used_filteres{$filter} = 1;
				$duplicate++;
			}					
		}
		if ($duplicate > 1)
		{
			return 0;
		}
		
	}

	my @correct_genotypes = sort( keys %filters_results );
	if (! @correct_genotypes)
	{
		return 0;
	}
	
	if ($num_filters == 3)
	{
		if ( scalar(@correct_genotypes) == $num_filters)
		{
			return 1;
		}
		else
		{
			return 0;
		}
	}
	elsif ( $num_filters == 2 )
	{
		if ( $oneFiltGroup_twoGenotypes == 1  )
		{
			if ( $number_of_genotypes > 2 and scalar( @correct_genotypes ) == 1 )
			{
				if ( $correct_genotypes[0] eq "[-/-]" or $correct_genotypes[0] eq "[+/+]" )
				{
					my $total_match = 0;
					my $total_genotypes = 0;
					my $selected_filter;
					foreach my $genot ( grep { $_ =~ /^\[/ } keys %{$$data_filtered{genot}} )
					{
						next if ( $genot eq $correct_genotypes[0] );

						foreach my $filter ( keys %$cfw_groups_thread )
						{
							if ($filters_results{$genot})
							{
								next if ( $filter == $filters_results{$genot} );
							}
							
							$selected_filter = $filter;
							$total_genotypes += scalar( split("\t", $$data_filtered{genot}{$genot}) );
							if ( $filters_matched_indv{$genot}{$filter} )
							{
								$total_match += $filters_matched_indv{$genot}{$filter};
							}										
						}									
					}
					
					my $err_percent = 1 - ( $total_match / scalar( @{$$cfw_groups_thread{$selected_filter}{indv}} ) );
					
					if ( $err_percent <= $$cfw_groups_thread{$selected_filter}{max_err}  )
					{
						return 1;
					}								
				}
				else
				{
					return 0;
				}
			}
			elsif ( scalar( @correct_genotypes ) == 2 )
			{
				return 1;
			}
			else
			{
				return 0;
			}						
		}
		else
		{
			if ( scalar( @correct_genotypes ) == 2)
			{
				return 1;
			}
			else
			{
				return 0;
			}
		}
	}				
}



#-----------------------------------------------------------------------------#
# The subroutine writes CAPS markers that have passed the genotype filtration #
#-----------------------------------------------------------------------------#
sub write2file_filtered_CAPS
{
	my $data_filtered = $_[0];
	
	my $Ofh;
	if ( !open $Ofh, '>>', $cfw_gf_output_file )
	{
		return ($!, 6);
	}
	
	print $Ofh "$$data_filtered{id}\n";
	print $Ofh "$$data_filtered{enz}\n";
	print $Ofh "$$data_filtered{ref}\n";
	
	foreach my $alt ( split(":::", $$data_filtered{alt}) )
	{
		print $Ofh "$alt\n";
	}
	
	foreach my $genot ( split(":::", $$data_filtered{genot_print}) )
	{
		print $Ofh "$genot\n";
	}
	
	close $Ofh;
	
	return ("", 1);
}


#--------------------------------------------------------#
# The subroutine converts the VCF file into the v2c file #
#--------------------------------------------------------#
sub vcf2snps
{
	my $genotypesLength = $_[0];
	my $reference_file_name_tmp = $_[1];
	my $error = 0;
	my $v2c_file_name_tmp = $raw_vcf_file_name;
	my %chrom_offset = (); # Stores data from the reference index file
	my $fh;
	
	$line_vcf = 0;
	@sequencesNotPresentInRef = ();
	$snpsNo = 0;
	@markersOnTheEdge = (); 
	
	$v2c_file_name_tmp =~ s/.*[\\\/]//g;
	$v2c_file_name_tmp =~ s/\..+$/\.v2c/;
	
	$reference_file_name_tmp =~ s/\..+$//;
	$reference_file_name_tmp =~ s/.*[\\\/]//g;
	
	if ( !open $fh, '<', $working_dir . $reference_file_name_tmp . ".index" )
	{
		return($!,3);
	}
	
		while (<$fh>)
		{
			chomp $_;			
			my @data = split("\t", $_);
			if ($data[0] =~ /^#/) { next }
			
			$data[0] =~ s/>//;
			$data[0] =~ s/ .+//g;
			$chrom_offset{$data[0]} = [($data[1],$data[2],$data[3],$data[4])];	
		}
	close $fh;
	
	my $chrom_offset = \%chrom_offset;
	
	
	if (-e $working_dir . $v2c_file_name_tmp) {unlink $working_dir . $v2c_file_name_tmp}
	
	if ( !open $fh, '<', $raw_vcf_file_name )
	{
		return($!,4);
	}
	
	my $Ofh;
	if ( !open $Ofh, '>>', $working_dir . $v2c_file_name_tmp )
	{
		return($!,5);
	}

	print $Ofh "#$reference_md5\n";
	print $Ofh "#$vcf_md5\n";
	
	while (my $bier = <$fh>)
	{
		chomp $bier;
		my @input = (split("\t", $bier));
		
		if ($input[0] =~ /^#/) { next }
		
		$line_vcf++;
		my $RefSnpLen = scalar(split("", $input[3])); # Reference variant length
		my $ref_snp = $input[3]; # Reference variant sequence
		my $AltSnp = $input[4];
		my @AltSnp = split(",", $AltSnp); # Alternative variant sequence(s)
		my $from = $input[1] - $genotypesLength;
		my $snp_min = $input[1] - 1; # The number of the nucleotide located just before the variant
		my $to = $input[1] + $genotypesLength + $RefSnpLen - 1;
		my $RefSnp_plus = $input[1] + $RefSnpLen; # The number of the nucleotide located just behind the variant
		
		if (!$chrom_offset->{$input[0]})
		{
			push @sequencesNotPresentInRef, ">" . $input[0];
			next;
		}
		elsif ($from < 0 or $to > $chrom_offset->{$input[0]}[3])
		{
			push @markersOnTheEdge, ">" . $input[0] . ":" . $input[1];
			next;
		}
				
		foreach my $snp (@AltSnp)
		{
			$snp =~ s/\s//g;
		}		
		
		my $AltSnpNo = @AltSnp;
		
		my @AltSnpLen = ();
		for (my $i = 0; $i < $AltSnpNo; $i++)
		{
			push @AltSnpLen, scalar(split("", $AltSnp[$i]));
		}

		my @genotypes = ();
		
		for (my $i = 9; $i < scalar(@input); $i++)
		{
			my @genotype = split(":", $input[$i]);
			my @alleles = split(/[\|\/]/, $genotype[0]);
			
			if ( $alleles[0] eq "." or $alleles[1] eq "." )
			{
				push @genotypes, $genotype[0];
			}
			else
			{
				if ( $alleles[0] <= $alleles[1] )
				{
					push @genotypes, $genotype[0];
				}
				else
				{
					push @genotypes, $alleles[1] . "/" . $alleles[0];
				}
			}			
		}

		my @seqExtractorOut = seqExtractor($reference_file_name,"$input[0]:$from-$snp_min",$chrom_offset);		
		
		print $Ofh ">$input[0]:$input[1]\n";
		print $Ofh join("\t",@genotypes)."\n";
		
		my $extracted_sequence = $seqExtractorOut[0];	

		print $Ofh "$ref_snp,$RefSnpLen,$extracted_sequence$ref_snp";
		
		@seqExtractorOut = seqExtractor($reference_file_name,"$input[0]:$RefSnp_plus-$to",$chrom_offset);
		
		$extracted_sequence = $seqExtractorOut[0];

		print $Ofh "$extracted_sequence";

		for (my $i = 0; $i < $AltSnpNo; $i++)
		{
			my $AltSnp_plus = $input[1] + $AltSnpLen[$i];
			$to = $input[1] + $genotypesLength + $AltSnpLen[$i] - 1;
			
			@seqExtractorOut = seqExtractor($reference_file_name,"$input[0]:$from-$snp_min",$chrom_offset);
			
			$extracted_sequence = $seqExtractorOut[0];

			print $Ofh "\t$AltSnp[$i],$AltSnpLen[$i],$extracted_sequence$AltSnp[$i]";
			
			@seqExtractorOut = seqExtractor($reference_file_name,"$input[0]:$AltSnp_plus-$to",$chrom_offset);
		
			$extracted_sequence = $seqExtractorOut[0];

			print $Ofh "$extracted_sequence";
		}
		print $Ofh "\n";
		
		$snpsNo++;
	}
	return ("",1);
}


#------------------------------------------------------------#
# The subroutine finds restriction sites within the v2c file #
#------------------------------------------------------------#
sub caps_miner
{
	my $seq_len = $_[0];
	my $reference_file_name_tmp = $_[1];
	$reference_file_name_tmp =~ s/\..+$//;
	$reference_file_name_tmp =~ s/.*[\\\/]//g;
	my $seqName;
	my @genotypes;
	my %genotypes;
	my @seq;
	my $linie_ile = scalar(@linie);
	my %chrom_offset = ();
	my $fh;
	
	$enzyme_zip = 0;
	$stop = 0;
	@markersOnTheEdge = ();	
	$numberOfSNPs = 0;

	if (-e "$working_dir" . "caps_markers.txt") { unlink "$working_dir" . "caps_markers.txt" }

	if ( !open $fh, '<', $working_dir . $reference_file_name_tmp . ".index" )
	{
		return($!, 3);
	}
		while (<$fh>)
		{
			chomp $_;
			my @data = split("\t", $_);
			if ($data[0] =~ /^#/) { next }
			$data[0] =~ s/>//;
			$data[0] =~ s/ .+//g;
			$chrom_offset{$data[0]} = [($data[1],$data[2],$data[3],$data[4])];
		}
	close $fh;
	my $chrom_offset = \%chrom_offset;
	
	my @selected_enz_names_tmp;
	my %selected_enz_names_mix;
	foreach my $enzyme_mix (keys %enzymes_for_analysis)
	{
		L: foreach my $enzyme (@selected_enz_names)
		{	
			foreach my $single_enz ( split(",", $enzyme_mix) )
			{
				if ($single_enz eq $enzyme)
				{
					push @selected_enz_names_tmp, $enzyme;
					$selected_enz_names_mix{$enzyme} = $enzyme_mix;
					last L;
				}
			}					
		}
		$enzyme_zip++;
	}
	@selected_enz_names_tmp = uniq(@selected_enz_names_tmp);
	
	$actualSNPNo = 0;	
	
	if ( !open $fh, '>>', "$working_dir" . "caps_markers.txt" )
	{
		return($!, 4);
	}
	print $fh "#VCF2CAPS_output_file\n\n";
	
	my $SNP;
	if ( !open $SNP, '<', $v2c_file_name )
	{
		return($!, 5);
	}
	
	while (my $bierSNP = <$SNP>)
	{
		if ($stop == 1) { return ($numberOfSNPs, 2) }
		
		chomp $bierSNP;
		
		my $IndvSeq = [];
		my @IndvSeq_snpLen_snp = ();
		my @IndvSeq_snp = ();
		
		if ($bierSNP =~ /^#/) { next }
		elsif ($bierSNP =~ /^>/)
		{
			$seqName = $bierSNP;
			$actualSNPNo++;
		}
		elsif ($bierSNP =~ /^[0-9\.]/)
		{
			@genotypes = split("\t",$bierSNP);
			
			for (my $i = 0; $i < scalar(@linie); $i++)
			{
				$genotypes{$linie[$i]} = $genotypes[$i];
			}
		}
		else
		{
			@seq = split("\t",$bierSNP);
			my $SNP_alleles_No = @seq;

			my $enz_No = @selected_enz_names_tmp;

			for (my $i = 0; $i < $SNP_alleles_No; $i++)
			{
				my @input = (split(",",$seq[$i]));
				push @IndvSeq_snpLen_snp, $input[1];
				push @IndvSeq_snp, $input[0];
				my @IndvSeq_tmp = split("",$input[2]);
				push @$IndvSeq, [@IndvSeq_tmp];
			}			
			
			for (my $i = 0; $i < $enz_No; $i++)
			{			
				my $enzLength = split("", ( split("\t",$enzymes_db{$selected_enz_names_tmp[$i]}) )[1] );

				$enzLength = $enzLength - 0;

				my $cc = 0;
				my @partSeq = ([],[]);
				my @partSeq_size = ();
				my $numberOfMatches = 0;
				my $numberOfMatches_all = 0;
				my %matches1;
				my @partSeq1 = "";
				
				my ($to, $count) = ($enzLength,0);
				
				for (my $seqID = 0; $seqID < $SNP_alleles_No; $seqID++)
				{
					$cc = 0;
					my $IndvSeq_Len = (scalar(@{$IndvSeq->[$seqID]}));
					my $from = int( ( ( $IndvSeq_Len - $IndvSeq_snpLen_snp[$seqID] ) + 1 ) / 2 ) - $enzLength;

					for (my $c = 1; $c < $enzLength + ($IndvSeq_snpLen_snp[$seqID]); $c++)
					{
						my @part = ();
					
						for (my $ii = $from + $c; $ii < $from + $enzLength + $c; $ii++)
						{

							my $input = ${@$IndvSeq[$seqID]}[$ii];
							push @part, $input;
							
						}
						
						$partSeq[$seqID][$cc] = join("",@part);
						$cc++;											
					}
				}

				for (my $i = 0; $i < $SNP_alleles_No; $i++)
				{
					my $input = @{$partSeq[$i]};
					@partSeq_size = (@partSeq_size, $input);
				}

				my $enzyme_recogn_seq = ( split("\t",$enzymes_db{$selected_enz_names_tmp[$i]}) )[1];
				$enzyme_recogn_seq = uc $enzyme_recogn_seq;
				
				my $enzREGEX = enzREGEX( $enzyme_recogn_seq );
				
				my $regex_inv = regex_inv($enzREGEX);
				if ($regex_inv eq $enzREGEX)
				{
					for (my $c = 0; $c < $SNP_alleles_No; $c++)
					{ 
						for (my $i = 0; $i < $partSeq_size[$c]; $i++)
						{	
							my $partSeq_upper = uc $partSeq[$c][$i];
							if ($partSeq_upper =~ /$enzREGEX/)
							{								
								$numberOfMatches++;
								$numberOfMatches_all++;
							}
						}
						$matches1{$c} = $numberOfMatches;
						$numberOfMatches = 0;
					}
				}
				else
				{
					for (my $c = 0; $c < $SNP_alleles_No; $c++)
					{ 
						for (my $i = 0; $i < $partSeq_size[$c]; $i++)
						{	
							my $partSeq_upper = uc $partSeq[$c][$i];
							if ($partSeq_upper =~ /$enzREGEX/ or $partSeq_upper =~ /$regex_inv/)
							{								
								$numberOfMatches++;
								$numberOfMatches_all++;
							}
						}
						$matches1{$c} = $numberOfMatches;
						$numberOfMatches = 0;
					}
				}

				my $DigestedGenotypes_No_homo = 0;
				my $DigestedGenotypes_No_het = 0;
				my $nullGenotypes = 0;

				for (my $i = 0; $i < $linie_ile; $i++)
				{
					my $genotypp = $genotypes{$linie[$i]};

					my @genotypp_indv = split("/", $genotypp);
					
					if ($genotypp_indv[0] ne "\.")
					{
						if ($genotypp_indv[0] == $genotypp_indv[1])
						{
							if ($matches1{$genotypp_indv[0]} > 0)
							{
								$DigestedGenotypes_No_homo++;
							}
						}
						elsif ($genotypp_indv[0] != $genotypp_indv[1])
						{
							my $check_digest = 0;
							my $check_NO_digest = 0;
							foreach my $genot (@genotypp_indv)
							{
								if ($matches1{$genot} > 0)
								{
									$check_digest++;
								}
								else
								{
									$check_NO_digest++;
								}
							}
							
							if ($check_digest > 0 and $check_NO_digest == 0)
							{
								$DigestedGenotypes_No_homo++;
							}
							elsif ($check_digest > 0 and $check_NO_digest > 0)
							{
								$DigestedGenotypes_No_het++;
							}
						}
					}
					else
					{
						$nullGenotypes++;
					}
				}
				
				my @genotypes_uniq = uniq(@genotypes);
				my @genotypes_uniq_tmp = grep { $_ ne "\.\/\." } @genotypes_uniq;
				
				if ($numberOfMatches_all > 0 and ( $DigestedGenotypes_No_homo > 0 or $DigestedGenotypes_No_het > 0 ) and $DigestedGenotypes_No_homo < ( $linie_ile - $nullGenotypes ) and scalar(@genotypes_uniq_tmp) > 1 )
				{ 
					$numberOfSNPs++;

					my ($chrom,$snp) = split(":",$seqName);
					my $snp_min = $snp - 1;
					$chrom =~ s/>//;
					my $taker_from = $snp - $seq_len;
					my $snp_plus = $snp + $IndvSeq_snpLen_snp[0];
					my $taker_to = $snp + $seq_len + $IndvSeq_snpLen_snp[0] - 1;
					
					if ( $taker_from < 0 or $taker_to > $chrom_offset->{$chrom}[3] )
					{
						push @markersOnTheEdge, $seqName;
						$numberOfSNPs--;
						next;
					}

					print $fh "$seqName\n";
					my $enz_seq = ( split("\t",$enzymes_db{$selected_enz_names_tmp[$i]}) )[-1];
					my $tmp = $selected_enz_names_mix{$selected_enz_names_tmp[$i]};
					print $fh "$tmp\t$enz_seq\n";
					
					for (my $i = 0; $i < $SNP_alleles_No; $i++)
					{
						if ($i == 0)
						{
							print $fh "ref\t[$i]\t$IndvSeq_snp[$i]\t";
						}
						else
						{
							print $fh "alt\t[$i]\t$IndvSeq_snp[$i]\t";
						}
						
						if ($matches1{$i} == 0)
						{
							print $fh "[-]\t";
						}
						else
						{
							print $fh "[+]\t";
						}
						
						my @taker = seqExtractor($reference_file_name,"$chrom:$taker_from-$snp_min",$chrom_offset);
						$taker[0] = lc $taker[0];
						print $fh $taker[0] . $IndvSeq_snp[$i];
						
						@taker = seqExtractor($reference_file_name,"$chrom:$snp_plus-$taker_to",$chrom_offset);
						$taker[0] = lc $taker[0];
						print $fh "$taker[0]\n";
					}
					
					foreach my $genotype_uniq ( sort{$b cmp $a} @genotypes_uniq )
					{
						print $fh "$genotype_uniq\t";
						print $fh "[";
						
						my $c = 0;
						foreach my $indv_allele ( split("/", $genotype_uniq) )
						{
							if ($indv_allele eq "\.")
							{
								print $fh "?";
							}
							else
							{
								if ($matches1{$indv_allele} == 0)
								{
									print $fh "-";
								}
								else
								{
									print $fh "+";
								}
							}
							
							if ($c == 0)
							{
								print $fh "/";
							}
							else
							{
								print $fh "]";
								
								foreach my $indv_name (sort{$a cmp $b} keys %genotypes)
								{										
									if ($genotypes{$indv_name} eq $genotype_uniq)
									{
										print $fh "\t$indv_name";
									}
								}
							}
							
							$c++;
						}
						print $fh "\n";
					}
					print $fh "\n";					
				}
			}
		}
	}
	
	close $fh;
	return ($numberOfSNPs,1);
}

sub singleCutSite
{
	my $EnzSeq = "";
	my @EnzSeqSingleNucl;
	my $EnzSeqSingleNuclSize;
	my %AllSeq;
	my %AllSeq_NoOfRecognSeq;
	my $preAllSeq;
	my @input;
	my @output = ();
	my $fh;
	my $Ofh;
	$numberOfSNPsBefore = 0;
	$numberOfSNPsAfter = 0;
	$cfw_scf_output_file = $cfw_scf_input_file;
	$cfw_scf_output_file =~ s/\.txt$/_scf\.txt/;
	$cfw_scf_output_file =~ s/.*[\\\/]//g;
	$cfw_scf_output_file = $working_dir . $cfw_scf_output_file;
	$stop = 0;
	
	if (-e $cfw_scf_output_file) { unlink $cfw_scf_output_file }

	if ( !open $Ofh, '>>', $cfw_scf_output_file )
	{
		return($!, 3);
	}
	
	print $Ofh "#VCF2CAPS_output_file\n\n";
	
	if ( !open $fh, '<', $cfw_scf_input_file )
	{
		return($!, 4);
	}

		while (my $bier = <$fh>)
		{
			if ($stop == 1) { return ("",2) }
			
			chomp $bier;
			if ($bier =~ /^>/)
			{
				if ( scalar(keys %AllSeq) != 0 )
				{
					my $check = 0;
					foreach my $snp_allele_ID (keys %AllSeq)
					{
						if ($AllSeq_NoOfRecognSeq{$snp_allele_ID} > 1)
						{
							$check++;
						}
					}
					
					if ($check == 0)
					{
						foreach my $line (@output)
						{
							print $Ofh "$line\n";
						}
						
						print $Ofh "\n";
						$numberOfSNPsAfter++;
					}
				}
								
				%AllSeq = ();
				%AllSeq_NoOfRecognSeq = ();
				@output = ();
				
				push @output, $bier;
				$numberOfSNPsBefore++;
			}
			elsif ($bier =~ /^[A-Z]/ and scalar(split("\t", $bier)) == 2 )
			{
				@input = split("\t",$bier);
				$EnzSeq = $input[1];
				$EnzSeq =~ s/[_']//g;
				$EnzSeq =~ s/(^[nN]+)|([nN]+$)//g;
				push @output, $bier;
			}
			elsif ($bier =~ /^[ra]/ and scalar(split("\t", $bier)) == 5)
			{
				push @output, $bier;
				@input = split("\t",$bier);
				$preAllSeq = $input[4];
				$preAllSeq = uc $preAllSeq;
				my $snp_ID = $input[1];
				$snp_ID =~ s/[\[\]]//g;
				$AllSeq{$snp_ID} = $preAllSeq;
			}
			elsif ($bier =~ /^[0-9\.]/ and scalar(split("\t", $bier)) >= 3)
			{
				push @output, $bier;
			}
			
			$EnzSeq = uc $EnzSeq;
			my $enzREGEX = enzREGEX($EnzSeq);
			
			@EnzSeqSingleNucl = split("",$EnzSeq);
			$EnzSeqSingleNuclSize = @EnzSeqSingleNucl; 
			
			foreach my $snp_allele_ID (keys %AllSeq)
			{
				my $NoOfRecognSeq = 0;
				my @Seq_singleNucl = split("", $AllSeq{$snp_allele_ID});
				my $Seq_singleNucl_size = @Seq_singleNucl;
				
				my $regex_inv = regex_inv($enzREGEX);
				if ($regex_inv eq $enzREGEX)
				{
					for (my $i = 0; $i < $Seq_singleNucl_size - $EnzSeqSingleNuclSize; $i++)
					{
						my @seqPart = @Seq_singleNucl[$i..( $i + ($EnzSeqSingleNuclSize - 1) )];
						my $seqPart = join("",@seqPart);
						my $seqPart_upper = uc $seqPart;
						if (uc $seqPart_upper =~ /$enzREGEX/)
						{
							$NoOfRecognSeq++;
						}
					}
				}
				else
				{
					for (my $i = 0; $i < $Seq_singleNucl_size - $EnzSeqSingleNuclSize; $i++)
					{
						my @seqPart = @Seq_singleNucl[$i..( $i + ($EnzSeqSingleNuclSize - 1) )];
						my $seqPart = join("",@seqPart);
						my $seqPart_upper = uc $seqPart;
						if (uc $seqPart_upper =~ /$enzREGEX/ or uc $seqPart_upper =~ /$regex_inv/)
						{
							$NoOfRecognSeq++;
						}
					}
				}

				$AllSeq_NoOfRecognSeq{$snp_allele_ID} = $NoOfRecognSeq;
			}
		}
	close $fh;	
	
	my $check = 0;
	foreach my $snp_allele_ID (keys %AllSeq)
	{
		if ($AllSeq_NoOfRecognSeq{$snp_allele_ID} > 1)
		{
			$check++;
		}
	}
	
	if ($check == 0)
	{
		foreach my $line (@output)
		{
			print $Ofh "$line\n";
		}
		
		print $Ofh "\n";
		$numberOfSNPsAfter++;
	}
	
	close $Ofh;
	return ("",1);
}

sub enzREGEX
{
	my $regex = "";
	my ($Seq) = @_;
	my @SeqIndv = split("",$Seq);
	my $SeqIndvL = @SeqIndv;
	for (my $i=0;$i<$SeqIndvL;$i++) {
		if ($SeqIndv[$i] =~ /[ATGC]/) {$regex = $regex.$SeqIndv[$i]}
		elsif ($SeqIndv[$i] eq "R") {$regex = $regex."[AG]"}
		elsif ($SeqIndv[$i] eq "Y") {$regex = $regex."[CT]"}
		elsif ($SeqIndv[$i] eq "S") {$regex = $regex."[CG]"}
		elsif ($SeqIndv[$i] eq "W") {$regex = $regex."[AT]"}
		elsif ($SeqIndv[$i] eq "K") {$regex = $regex."[GT]"}
		elsif ($SeqIndv[$i] eq "M") {$regex = $regex."[AC]"}
		elsif ($SeqIndv[$i] eq "B") {$regex = $regex."[CGT]"}
		elsif ($SeqIndv[$i] eq "D") {$regex = $regex."[AGT]"}
		elsif ($SeqIndv[$i] eq "H") {$regex = $regex."[ACT]"}
		elsif ($SeqIndv[$i] eq "V") {$regex = $regex."[ACG]"}
		elsif ($SeqIndv[$i] eq "N") {$regex = $regex."[ACGT]"}
	}
	return $regex;
}

sub uniq {
	my %seen;
	return grep {!$seen{$_}++} @_;
}

sub seqExtractor
{
	my ($fh,$Ofh);
	my $error_code = 1;
	my $input = $_[0];
	my $chrom_and_coord = $_[1];
	my $chrom = (split(":", $chrom_and_coord))[0];
	my $coord = (split(":", $chrom_and_coord))[1];
	my $from_nucl = (split("-", $coord))[0];
	my $to_nucl = (split("-", $coord))[1];
	my $chrom_offset_ = $_[2];
	
	if (!$chrom_offset_->{$chrom})
	{
		return ($chrom,3);
	}
	elsif ($from_nucl eq "" or $to_nucl > $chrom_offset_->{$chrom}[3])
	{
		return ($chrom,4);
	}
	
	open $fh, '<', $input or die "Cannot open the file $input\n";
	
		my $chrom_offset = $chrom_offset_->{$chrom}[0];
		my $line_len = $chrom_offset_->{$chrom}[1];
		my $last_line_len = $chrom_offset_->{$chrom}[2];
		my $chrom_len = $chrom_offset_->{$chrom}[3];

		if (!defined($chrom_offset))
		{
			$error_code = 3;
			return ($chrom,$error_code);
		}
		
		my $from_offset = 0;
		my $offset = 0;
		
		if ($from_nucl <= ($chrom_len - $last_line_len + 1) or $chrom_len == $last_line_len)
		{
			$from_offset = $chrom_offset + $from_nucl - 1 + int($from_nucl / $line_len);
			$offset = $to_nucl - $from_nucl + 1 + (int($to_nucl / $line_len) - int($from_nucl / $line_len));
			if (($from_nucl % $line_len) == 0) {$from_offset -= 1; $offset += 1}
		}
		else
		{
			$from_offset = $chrom_offset + $from_nucl + int($from_nucl / $line_len);
			$offset = $to_nucl - $from_nucl + (int($to_nucl / $line_len) - int($from_nucl / $line_len));
			if (($from_nucl % $line_len) != 0) {$from_offset -= 1; $offset += 1}
		}

		my $line = "";
		seek($fh, $from_offset, 0);
		read($fh, $line, $offset);
		$line =~ s/\n//g;

		return ($line,$error_code);

	close $fh;
}

sub regex_inv
{
	my $input = $_[0];
	my @data = split("", $input);
	my @data_inv = ();
	
	my @temp = ();
	my $check = 0;
	for (my $i = scalar(@data) - 1; $i >= 0; $i-- )
	{
		if ($data[$i] =~ /\]/)
		{
			$check = 1;
			next;
		}
		elsif ($data[$i] =~ /\[/)
		{
			my @temp_sorted = sort { $a cmp $b } @temp;
			push @temp_sorted, ']';
			unshift @temp_sorted, '[';
			push @data_inv, join("", @temp_sorted);
			
			@temp = ();
			$check = 0;
		}
		elsif ($check == 1 and $data[$i] =~ /[^\[\]]/ )
		{
			unshift @temp, invert($data[$i]);
		}
		else
		{
			push @data_inv, invert($data[$i]);
		}
	}
	
	return join("", @data_inv);
}

sub invert
{
	my $input = $_[0];
	my $output = "";
	
	if ($input eq 'A') { $output = 'T' }
	elsif ($input eq 'T') { $output = 'A' }
	elsif ($input eq 'C') { $output = 'G' }
	elsif ($input eq 'G') { $output = 'C' }
	else { $output = $input }
	
	return $output;
}

sub curr_time
{
	my ($sec, $min, $hour, $mday, $mon, $year, $wday, $yday, $isdst) = localtime(time);
	my $formated_time = sprintf("%02d:%02d:%02d", $hour, $min, $sec);
	$terminal->insert('end', "[" . $formated_time . "] ", 'date');
}

sub warning_dial
{
	my ($title, $message) = split(':::', $_[0]);
	
	if (defined $title and defined $message)
	{
		$mw->messageBox(
			-title => $title,
			-message => $message,
			-type => 'OK',
			-icon => 'error'
		);
	}
}
