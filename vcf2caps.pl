use strict;
use warnings;
use Tk;
use Tk::FBox;
use Tk::Animation;
use threads;
use threads::shared;
use Tk::Listbox;
use Tk::Pane;
use Tk::ProgressBar;
use LWP::Simple;
use Data::Dumper;


$|++;

my $enzyme_file_name:shared;
my $reference_file_name:shared;
my $vcf_file_name:shared;
my $raw_vcf_file_name:shared;
my $die:shared = 0;
my $die_confirm:shared = 0;
#my %enzymes_db:shared;
my $jobID:shared = 0;
my %enzymes_db:shared;
my $line_vcf:shared = 0;
my $line:shared = 0;
my @reference_analysis_results:shared = (0);
my @enzyme_analysis_results:shared = (0);
my %vcf_analysis_results:shared = (err_code => 0);
my %raw_vcf_analysis_results:shared = (err_code => 0);
my @allEnzymesNames:shared;
my @snp2caps_results:shared;
my @singleCutSite_results:shared;
my $numberOfSNPs:shared = 0;
my $fancy_SNPNo:shared = 0;
my $capsMining_percent;
my $singleCutSite_percent;
my @linie:shared;
my $numberOfSNPsAfter:shared;
my $numberOfSNPsBefore:shared = 0;
my %enzymes_for_analysis:shared;
my $enzyme_zip:shared = 0;
my $enzyme_zip_stop:shared = 0;
my $reference_start_analysis:shared = 0;
my $enzyme_start_analysis:shared = 0;
my $vcf_start_analysis:shared = 0;
my @sequencesNotPresentInRef:shared = ();
my $snpsNo:shared = 0;
my $stop:shared;
my @markersOnTheEdge:shared = ();
my $log_first_use = 0;
my $genome_exists = 0;

my $fh;
my $out;
my @selected_enz_names:shared = ();
my $selected_enz_name;
my @samtools_error_code:shared = ("",0);

my $enzyme_check;
my $enzyme_check_status;
my $reference_check;
my $reference_check_status;
my $vcf_check;
my $vcf_check_status;
my $raw_vcf_check;
my $raw_vcf_check_status;
my $snps_seq_len:shared = 40;
my $output_seq_len:shared = 500;
my $onlyCommercially = 0;
my $polymorphicSNPsOnly:shared = 1;
my $cutters4 = 1;
my $cutters5 = 1;
my $cutters6 = 0;
my $custom = 0;
my $custom_value = "";
my $iso_state = 0;
my $comp_state = 0;
my $working_dir:shared = "";

threads->create( \&work )->detach();

my $mw = MainWindow->new();
$mw->title("vcf2caps v2.0");
$mw->minsize( qw(500 300) );

my $L_center_col2_1_entry;
my $L_center_col2_2_entry;

##About START
$mw->fontCreate('title', -family => 'arial', -size => 12, -weight => 'bold');
$mw->fontCreate('text', -family => 'arial', -size => 8);
$mw->fontCreate('text_b', -family => 'arial', -size => 8, -weight => 'bold');
$mw->fontCreate('hyper', -family => 'arial', -size => 8);
my $about_window = $mw->Toplevel(-title => 'About vcf2caps');
$about_window->resizable(0,0);
my $about_text = $about_window->Text(-cursor => 'left_ptr', -width => 50, -height => 14, -insertwidth => 0, -insertontime => 0, -insertofftime => 0, -wrap => 'word', -background => 'gray95')->pack();
$about_text->tagConfigure('title_center', -justify => 'center', -font => 'title');
$about_text->tagConfigure('text_center', -justify => 'center', -font => 'text');
$about_text->tagConfigure('hyperlink', -underline => 0, -font => 'hyper', -foreground => 'blue', -justify => 'center');
$about_text->tagBind('hyperlink', '<Any-Enter>' => sub { $about_text->tagConfigure('hyperlink', -underline => 1, -font => 'hyper') } );
$about_text->tagBind('hyperlink', '<Any-Leave>' => sub { $about_text->tagConfigure('hyperlink', -underline => 0, -font => 'hyper') } );
$about_text->tagBind('hyperlink', '<Button-1>' => sub {
	open_hyperlink("http://google.pl");
} );
$about_text->insert('end',"\n");
$about_text->insert('end',"vcf2caps v2.0\n", 'title_center');
$about_text->insert('end',"Copyright \x{00A9} 2018\n\n", 'text_center');
$about_text->insert('end',"Free, open-source CAPS mining software from VCF files.\n", 'text_center');
$about_text->insert('end', "http://google.pl\n\n", 'hyperlink');
$about_text->insert('end',"Agricultural University of Krakow, Poland\n\n", 'text_center');
$about_text->insert('end',"Programmed by Wojciech Weso\x{0142}owski\n\n", 'text_center');
$about_text->insert('end'," ", 'text_center');

my $button_OK = $about_text->Button(-text => 'OK', -width => 10, -command => sub { $about_window->withdraw } );
$about_text->windowCreate('text_center.last', -window => $button_OK, -align => 'center');
$about_window->protocol('WM_DELETE_WINDOW' => sub { $about_window->withdraw } );
$about_window->withdraw;
##About END

##Licence START
my $licence_window = $mw->Toplevel(-title => 'About vcf2caps');
my $licence_text = $licence_window->Scrolled('Text', -scrollbars => 'e', -padx => 10, -pady => 10, -cursor => 'left_ptr', -width => 100, -height => 28, -insertwidth => 0, -insertontime => 0, -insertofftime => 0, -wrap => 'word', -background => 'gray95')->pack(-fill => 'both', -expand => 1);
$licence_text->tagConfigure('title_center', -justify => 'center', -font => 'title');
$licence_text->tagConfigure('title_left', -font => 'title');
$licence_text->tagConfigure('title_left_s', -font => 'text_b');
$licence_text->tagConfigure('text', -font => 'text');
$licence_text->tagConfigure('text_center', -font => 'text', -justify => 'center');
$licence_text->tagConfigure('hyperlink', -underline => 0, -font => 'hyper', -foreground => 'blue', -justify => 'center');
$licence_text->tagBind('hyperlink', '<Any-Enter>' => sub { $licence_text->tagConfigure('hyperlink', -underline => 1, -font => 'hyper') } );
$licence_text->tagBind('hyperlink', '<Any-Leave>' => sub { $licence_text->tagConfigure('hyperlink', -underline => 0, -font => 'hyper') } );
$licence_text->tagBind('hyperlink', '<Button-1>' => sub {
	open_hyperlink("http://google.pl");
} );
$licence_text->insert('end',"\n");
$licence_text->insert('end',"GNU GENERAL PUBLIC LICENSE\n", 'title_center');
$licence_text->insert('end',"Version 3, 29 June 2007\n\n", 'text_center');
$licence_text->insert('end',"Copyright \x{00A9} 2007 Free Software Foundation, Inc. <https://fsf.org/>\n", 'text');
$licence_text->insert('end',"Everyone is permitted to copy and distribute verbatim copies of this license document, but changing it is not allowed.\n\n", 'text');
$licence_text->insert('end',"Preamble\n\n", 'title_left');
$licence_text->insert('end',"The GNU General Public License is a free, copyleft license for software and other kinds of works.

The licenses for most software and other practical works are designed to take away your freedom to share and change the works. By contrast, the GNU General Public License is intended to guarantee your freedom to share and change all versions of a program--to make sure it remains free software for all its users. We, the Free Software Foundation, use the GNU General Public License for most of our software; it applies also to any other work released this way by its authors. You can apply it to your programs, too.

When we speak of free software, we are referring to freedom, not price. Our General Public Licenses are designed to make sure that you have the freedom to distribute copies of free software (and charge for them if you wish), that you receive source code or can get it if you want it, that you can change the software or use pieces of it in new free programs, and that you know you can do these things.

To protect your rights, we need to prevent others from denying you these rights or asking you to surrender the rights. Therefore, you have certain responsibilities if you distribute copies of the software, or if you modify it: responsibilities to respect the freedom of others.

For example, if you distribute copies of such a program, whether gratis or for a fee, you must pass on to the recipients the same freedoms that you received. You must make sure that they, too, receive or can get the source code. And you must show them these terms so they know their rights.

Developers that use the GNU GPL protect your rights with two steps: (1) assert copyright on the software, and (2) offer you this License giving you legal permission to copy, distribute and/or modify it.

For the developers' and authors' protection, the GPL clearly explains that there is no warranty for this free software. For both users' and authors' sake, the GPL requires that modified versions be marked as changed, so that their problems will not be attributed erroneously to authors of previous versions.

Some devices are designed to deny users access to install or run modified versions of the software inside them, although the manufacturer can do so. This is fundamentally incompatible with the aim of protecting users' freedom to change the software. The systematic pattern of such abuse occurs in the area of products for individuals to use, which is precisely where it is most unacceptable. Therefore, we have designed this version of the GPL to prohibit the practice for those products. If such problems arise substantially in other domains, we stand ready to extend this provision to those domains in future versions of the GPL, as needed to protect the freedom of users.

Finally, every program is threatened constantly by software patents. States should not allow patents to restrict development and use of software on general-purpose computers, but in those that do, we wish to avoid the special danger that patents applied to a free program could make it effectively proprietary. To prevent this, the GPL assures that patents cannot be used to render the program non-free.

The precise terms and conditions for copying, distribution and modification follow.\n\n", 'text');
$licence_text->insert('end',"TERMS AND CONDITIONS\n\n", 'title_left');
$licence_text->insert('end',"0. Definitions.\n\n", 'title_left_s');
$licence_text->insert('end',"\x{201F}This License\x{201D} refers to version 3 of the GNU General Public License.

\x{201F}Copyright\x{201D} also means copyright-like laws that apply to other kinds of works, such as semiconductor masks.

\x{201F}The Program\x{201D} refers to any copyrightable work licensed under this License. Each licensee is addressed as \x{201F}you\x{201D}. \x{201F}Licensees\x{201D} and \x{201F}recipients\x{201D} may be individuals or organizations.

To \x{201F}modify\x{201D} a work means to copy from or adapt all or part of the work in a fashion requiring copyright permission, other than the making of an exact copy. The resulting work is called a \x{201F}modified version\x{201D} of the earlier work or a work \x{201F}based on\x{201D} the earlier work.

A \x{201F}covered work\x{201D} means either the unmodified Program or a work based on the Program.

To \x{201F}propagate\x{201D} a work means to do anything with it that, without permission, would make you directly or secondarily liable for infringement under applicable copyright law, except executing it on a computer or modifying a private copy. Propagation includes copying, distribution (with or without modification), making available to the public, and in some countries other activities as well.

To \x{201F}convey\x{201D} a work means any kind of propagation that enables other parties to make or receive copies. Mere interaction with a user through a computer network, with no transfer of a copy, is not conveying.

An interactive user interface displays \x{201F}Appropriate Legal Notices\x{201D} to the extent that it includes a convenient and prominently visible feature that (1) displays an appropriate copyright notice, and (2) tells the user that there is no warranty for the work (except to the extent that warranties are provided), that licensees may convey the work under this License, and how to view a copy of this License. If the interface presents a list of user commands or options, such as a menu, a prominent item in the list meets this criterion.\n\n", 'text');
$licence_text->insert('end',"1. Source Code.\n\n", 'title_left_s');
$licence_text->insert('end',"The \x{201F}source code\x{201D} for a work means the preferred form of the work for making modifications to it. \x{201F}Object code\x{201D} means any non-source form of a work.

A \x{201F}Standard Interface\x{201D} means an interface that either is an official standard defined by a recognized standards body, or, in the case of interfaces specified for a particular programming language, one that is widely used among developers working in that language.

The \x{201F}System Libraries\x{201D} of an executable work include anything, other than the work as a whole, that (a) is included in the normal form of packaging a Major Component, but which is not part of that Major Component, and (b) serves only to enable use of the work with that Major Component, or to implement a Standard Interface for which an implementation is available to the public in source code form. A \x{201F}Major Component\x{201D}, in this context, means a major essential component (kernel, window system, and so on) of the specific operating system (if any) on which the executable work runs, or a compiler used to produce the work, or an object code interpreter used to run it.

The \x{201F}Corresponding Source\x{201D} for a work in object code form means all the source code needed to generate, install, and (for an executable work) run the object code and to modify the work, including scripts to control those activities. However, it does not include the work's System Libraries, or general-purpose tools or generally available free programs which are used unmodified in performing those activities but which are not part of the work. For example, Corresponding Source includes interface definition files associated with source files for the work, and the source code for shared libraries and dynamically linked subprograms that the work is specifically designed to require, such as by intimate data communication or control flow between those subprograms and other parts of the work.

The Corresponding Source need not include anything that users can regenerate automatically from other parts of the Corresponding Source.

The Corresponding Source for a work in source code form is that same work.
\n", 'text');
$licence_text->insert('end',"2. Basic Permissions.\n\n", 'title_left_s');
$licence_text->insert('end',"All rights granted under this License are granted for the term of copyright on the Program, and are irrevocable provided the stated conditions are met. This License explicitly affirms your unlimited permission to run the unmodified Program. The output from running a covered work is covered by this License only if the output, given its content, constitutes a covered work. This License acknowledges your rights of fair use or other equivalent, as provided by copyright law.

You may make, run and propagate covered works that you do not convey, without conditions so long as your license otherwise remains in force. You may convey covered works to others for the sole purpose of having them make modifications exclusively for you, or provide you with facilities for running those works, provided that you comply with the terms of this License in conveying all material for which you do not control copyright. Those thus making or running the covered works for you must do so exclusively on your behalf, under your direction and control, on terms that prohibit them from making any copies of your copyrighted material outside their relationship with you.

Conveying under any other circumstances is permitted solely under the conditions stated below. Sublicensing is not allowed; section 10 makes it unnecessary.
\n", 'text');
$licence_text->insert('end',"3. Protecting Users' Legal Rights From Anti-Circumvention Law.\n\n", 'title_left_s');
$licence_text->insert('end',"No covered work shall be deemed part of an effective technological measure under any applicable law fulfilling obligations under article 11 of the WIPO copyright treaty adopted on 20 December 1996, or similar laws prohibiting or restricting circumvention of such measures.

When you convey a covered work, you waive any legal power to forbid circumvention of technological measures to the extent such circumvention is effected by exercising rights under this License with respect to the covered work, and you disclaim any intention to limit operation or modification of the work as a means of enforcing, against the work's users, your or third parties' legal rights to forbid circumvention of technological measures.\n\n", 'text');
$licence_text->insert('end',"4. Conveying Verbatim Copies.\n\n", 'title_left_s');
$licence_text->insert('end',"You may convey verbatim copies of the Program's source code as you receive it, in any medium, provided that you conspicuously and appropriately publish on each copy an appropriate copyright notice; keep intact all notices stating that this License and any non-permissive terms added in accord with section 7 apply to the code; keep intact all notices of the absence of any warranty; and give all recipients a copy of this License along with the Program.

You may charge any price or no price for each copy that you convey, and you may offer support or warranty protection for a fee.\n\n", 'text');
$licence_text->insert('end',"5. Conveying Modified Source Versions.\n\n", 'title_left_s');
$licence_text->insert('end',"You may convey a work based on the Program, or the modifications to produce it from the Program, in the form of source code under the terms of section 4, provided that you also meet all of these conditions:

  a) The work must carry prominent notices stating that you modified it, and giving a relevant date.
  
  b) The work must carry prominent notices stating that it is released under this License and any conditions added under section 7. This requirement modifies the requirement in section 4 to “keep intact all notices”.
  
  c) You must license the entire work, as a whole, under this License to anyone who comes into possession of a copy. This License will therefore apply, along with any applicable section 7 additional terms, to the whole of the work, and all its parts, regardless of how they are packaged. This License gives no permission to license the work in any other way, but it does not invalidate such permission if you have separately received it.
  
  d) If the work has interactive user interfaces, each must display Appropriate Legal Notices; however, if the Program has interactive interfaces that do not display Appropriate Legal Notices, your work need not make them do so.
  
A compilation of a covered work with other separate and independent works, which are not by their nature extensions of the covered work, and which are not combined with it such as to form a larger program, in or on a volume of a storage or distribution medium, is called an “aggregate” if the compilation and its resulting copyright are not used to limit the access or legal rights of the compilation's users beyond what the individual works permit. Inclusion of a covered work in an aggregate does not cause this License to apply to the other parts of the aggregate.\n\n", 'text');
$licence_text->insert('end',"6. Conveying Non-Source Forms.\n\n", 'title_left_s');
$licence_text->insert('end',"You may convey a covered work in object code form under the terms of sections 4 and 5, provided that you also convey the machine-readable Corresponding Source under the terms of this License, in one of these ways:

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
$licence_text->insert('end',"7. Additional Terms.\n\n", 'title_left_s');
$licence_text->insert('end',"\x{201F}Additional permissions\x{201D} are terms that supplement the terms of this License by making exceptions from one or more of its conditions. Additional permissions that are applicable to the entire Program shall be treated as though they were included in this License, to the extent that they are valid under applicable law. If additional permissions apply only to part of the Program, that part may be used separately under those permissions, but the entire Program remains governed by this License without regard to the additional permissions.

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
$licence_text->insert('end',"8. Termination.\n\n", 'title_left_s');
$licence_text->insert('end',"You may not propagate or modify a covered work except as expressly provided under this License. Any attempt otherwise to propagate or modify it is void, and will automatically terminate your rights under this License (including any patent licenses granted under the third paragraph of section 11).

However, if you cease all violation of this License, then your license from a particular copyright holder is reinstated (a) provisionally, unless and until the copyright holder explicitly and finally terminates your license, and (b) permanently, if the copyright holder fails to notify you of the violation by some reasonable means prior to 60 days after the cessation.

Moreover, your license from a particular copyright holder is reinstated permanently if the copyright holder notifies you of the violation by some reasonable means, this is the first time you have received notice of violation of this License (for any work) from that copyright holder, and you cure the violation prior to 30 days after your receipt of the notice.

Termination of your rights under this section does not terminate the licenses of parties who have received copies or rights from you under this License. If your rights have been terminated and not permanently reinstated, you do not qualify to receive new licenses for the same material under section 10.\n\n", 'text');
$licence_text->insert('end',"9. Acceptance Not Required for Having Copies.\n\n", 'title_left_s');
$licence_text->insert('end',"You are not required to accept this License in order to receive or run a copy of the Program. Ancillary propagation of a covered work occurring solely as a consequence of using peer-to-peer transmission to receive a copy likewise does not require acceptance. However, nothing other than this License grants you permission to propagate or modify any covered work. These actions infringe copyright if you do not accept this License. Therefore, by modifying or propagating a covered work, you indicate your acceptance of this License to do so.\n\n", 'text');
$licence_text->insert('end',"10. Automatic Licensing of Downstream Recipients.\n\n", 'title_left_s');
$licence_text->insert('end',"Each time you convey a covered work, the recipient automatically receives a license from the original licensors, to run, modify and propagate that work, subject to this License. You are not responsible for enforcing compliance by third parties with this License.

An \x{201F}entity transaction\x{201D} is a transaction transferring control of an organization, or substantially all assets of one, or subdividing an organization, or merging organizations. If propagation of a covered work results from an entity transaction, each party to that transaction who receives a copy of the work also receives whatever licenses to the work the party's predecessor in interest had or could give under the previous paragraph, plus a right to possession of the Corresponding Source of the work from the predecessor in interest, if the predecessor has it or can get it with reasonable efforts.

You may not impose any further restrictions on the exercise of the rights granted or affirmed under this License. For example, you may not impose a license fee, royalty, or other charge for exercise of rights granted under this License, and you may not initiate litigation (including a cross-claim or counterclaim in a lawsuit) alleging that any patent claim is infringed by making, using, selling, offering for sale, or importing the Program or any portion of it.\n\n", 'text');
$licence_text->insert('end',"11. Patents.\n\n", 'title_left_s');
$licence_text->insert('end',"A \x{201F}contributor\x{201D} is a copyright holder who authorizes use under this License of the Program or a work on which the Program is based. The work thus licensed is called the contributor's \x{201F}contributor version\x{201D}.

A contributor's \x{201F}essential patent claims\x{201D} are all patent claims owned or controlled by the contributor, whether already acquired or hereafter acquired, that would be infringed by some manner, permitted by this License, of making, using, or selling its contributor version, but do not include claims that would be infringed only as a consequence of further modification of the contributor version. For purposes of this definition, \x{201F}control\x{201D} includes the right to grant patent sublicenses in a manner consistent with the requirements of this License.

Each contributor grants you a non-exclusive, worldwide, royalty-free patent license under the contributor's essential patent claims, to make, use, sell, offer for sale, import and otherwise run, modify and propagate the contents of its contributor version.

In the following three paragraphs, a \x{201F}patent license\x{201D} is any express agreement or commitment, however denominated, not to enforce a patent (such as an express permission to practice a patent or covenant not to sue for patent infringement). To \x{201F}grant\x{201D} such a patent license to a party means to make such an agreement or commitment not to enforce a patent against the party.

If you convey a covered work, knowingly relying on a patent license, and the Corresponding Source of the work is not available for anyone to copy, free of charge and under the terms of this License, through a publicly available network server or other readily accessible means, then you must either (1) cause the Corresponding Source to be so available, or (2) arrange to deprive yourself of the benefit of the patent license for this particular work, or (3) arrange, in a manner consistent with the requirements of this License, to extend the patent license to downstream recipients. \x{201F}Knowingly relying\x{201D} means you have actual knowledge that, but for the patent license, your conveying the covered work in a country, or your recipient's use of the covered work in a country, would infringe one or more identifiable patents in that country that you have reason to believe are valid.

If, pursuant to or in connection with a single transaction or arrangement, you convey, or propagate by procuring conveyance of, a covered work, and grant a patent license to some of the parties receiving the covered work authorizing them to use, propagate, modify or convey a specific copy of the covered work, then the patent license you grant is automatically extended to all recipients of the covered work and works based on it.

A patent license is \x{201F}discriminatory\x{201D} if it does not include within the scope of its coverage, prohibits the exercise of, or is conditioned on the non-exercise of one or more of the rights that are specifically granted under this License. You may not convey a covered work if you are a party to an arrangement with a third party that is in the business of distributing software, under which you make payment to the third party based on the extent of your activity of conveying the work, and under which the third party grants, to any of the parties who would receive the covered work from you, a discriminatory patent license (a) in connection with copies of the covered work conveyed by you (or copies made from those copies), or (b) primarily for and in connection with specific products or compilations that contain the covered work, unless you entered into that arrangement, or that patent license was granted, prior to 28 March 2007.

Nothing in this License shall be construed as excluding or limiting any implied license or other defenses to infringement that may otherwise be available to you under applicable patent law.\n\n", 'text');
$licence_text->insert('end',"12. No Surrender of Others' Freedom.\n\n", 'title_left_s');
$licence_text->insert('end',"If conditions are imposed on you (whether by court order, agreement or otherwise) that contradict the conditions of this License, they do not excuse you from the conditions of this License. If you cannot convey a covered work so as to satisfy simultaneously your obligations under this License and any other pertinent obligations, then as a consequence you may not convey it at all. For example, if you agree to terms that obligate you to collect a royalty for further conveying from those to whom you convey the Program, the only way you could satisfy both those terms and this License would be to refrain entirely from conveying the Program.\n\n", 'text');
$licence_text->insert('end',"13. Use with the GNU Affero General Public License.\n\n", 'title_left_s');
$licence_text->insert('end',"Notwithstanding any other provision of this License, you have permission to link or combine any covered work with a work licensed under version 3 of the GNU Affero General Public License into a single combined work, and to convey the resulting work. The terms of this License will continue to apply to the part which is the covered work, but the special requirements of the GNU Affero General Public License, section 13, concerning interaction through a network will apply to the combination as such.\n\n", 'text');
$licence_text->insert('end',"14. Revised Versions of this License.\n\n", 'title_left_s');
$licence_text->insert('end',"The Free Software Foundation may publish revised and/or new versions of the GNU General Public License from time to time. Such new versions will be similar in spirit to the present version, but may differ in detail to address new problems or concerns.

Each version is given a distinguishing version number. If the Program specifies that a certain numbered version of the GNU General Public License “or any later version” applies to it, you have the option of following the terms and conditions either of that numbered version or of any later version published by the Free Software Foundation. If the Program does not specify a version number of the GNU General Public License, you may choose any version ever published by the Free Software Foundation.

If the Program specifies that a proxy can decide which future versions of the GNU General Public License can be used, that proxy's public statement of acceptance of a version permanently authorizes you to choose that version for the Program.

Later license versions may give you additional or different permissions. However, no additional obligations are imposed on any author or copyright holder as a result of your choosing to follow a later version.\n\n", 'text');
$licence_text->insert('end',"15. Disclaimer of Warranty.\n\n", 'title_left_s');
$licence_text->insert('end',"THERE IS NO WARRANTY FOR THE PROGRAM, TO THE EXTENT PERMITTED BY APPLICABLE LAW. EXCEPT WHEN OTHERWISE STATED IN WRITING THE COPYRIGHT HOLDERS AND/OR OTHER PARTIES PROVIDE THE PROGRAM “AS IS” WITHOUT WARRANTY OF ANY KIND, EITHER EXPRESSED OR IMPLIED, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE. THE ENTIRE RISK AS TO THE QUALITY AND PERFORMANCE OF THE PROGRAM IS WITH YOU. SHOULD THE PROGRAM PROVE DEFECTIVE, YOU ASSUME THE COST OF ALL NECESSARY SERVICING, REPAIR OR CORRECTION.\n\n", 'text');
$licence_text->insert('end',"16. Limitation of Liability.\n\n", 'title_left_s');
$licence_text->insert('end',"IN NO EVENT UNLESS REQUIRED BY APPLICABLE LAW OR AGREED TO IN WRITING WILL ANY COPYRIGHT HOLDER, OR ANY OTHER PARTY WHO MODIFIES AND/OR CONVEYS THE PROGRAM AS PERMITTED ABOVE, BE LIABLE TO YOU FOR DAMAGES, INCLUDING ANY GENERAL, SPECIAL, INCIDENTAL OR CONSEQUENTIAL DAMAGES ARISING OUT OF THE USE OR INABILITY TO USE THE PROGRAM (INCLUDING BUT NOT LIMITED TO LOSS OF DATA OR DATA BEING RENDERED INACCURATE OR LOSSES SUSTAINED BY YOU OR THIRD PARTIES OR A FAILURE OF THE PROGRAM TO OPERATE WITH ANY OTHER PROGRAMS), EVEN IF SUCH HOLDER OR OTHER PARTY HAS BEEN ADVISED OF THE POSSIBILITY OF SUCH DAMAGES.\n\n", 'text');
$licence_text->insert('end',"17. Interpretation of Sections 15 and 16.\n\n", 'title_left_s');
$licence_text->insert('end',"If the disclaimer of warranty and limitation of liability provided above cannot be given local legal effect according to their terms, reviewing courts shall apply local law that most closely approximates an absolute waiver of all civil liability in connection with the Program, unless a warranty or assumption of liability accompanies a copy of the Program in return for a fee.

END OF TERMS AND CONDITIONS\n\n", 'text');
$licence_text->insert('end',"How to Apply These Terms to Your New Programs\n\n", 'title_left');
$licence_text->insert('end',"If you develop a new program, and you want it to be of the greatest possible use to the public, the best way to achieve this is to make it free software which everyone can redistribute and change under these terms.

To do so, attach the following notices to the program. It is safest to attach them to the start of each source file to most effectively state the exclusion of warranty; and each file should have at least the “copyright” line and a pointer to where the full notice is found.

    \<one line to give the program's name and a brief idea of what it does.\>
    Copyright (C) \<year\>  <name of author>

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see \<https://www.gnu.org/licenses/\>.
Also add information on how to contact you by electronic and paper mail.

If the program does terminal interaction, make it output a short notice like this when it starts in an interactive mode:

    \<program\>  Copyright (C) \<year\>  <name of author\>
    This program comes with ABSOLUTELY NO WARRANTY; for details type 'show w'.
    This is free software, and you are welcome to redistribute it
    under certain conditions; type 'show c' for details.
The hypothetical commands 'show w' and 'show c' should show the appropriate parts of the General Public License. Of course, your program's commands might be different; for a GUI interface, you would use an \x{201F}about box\x{201D}.

You should also get your employer (if you work as a programmer) or school, if any, to sign a \x{201F}copyright disclaimer\x{201D} for the program, if necessary. For more information on this, and how to apply and follow the GNU GPL, see <https://www.gnu.org/licenses/\>.

The GNU General Public License does not permit incorporating your program into proprietary programs. If your program is a subroutine library, you may consider it more useful to permit linking proprietary applications with the library. If this is what you want to do, use the GNU Lesser General Public License instead of this License. But first, please read <https://www.gnu.org/licenses/why-not-lgpl.html\>.\n\n", 'text');


$licence_text->insert('end'," ", 'text_center');

my $button_OK_lic = $licence_text->Button(-text => 'Close', -width => 10, -command => sub { $licence_window->withdraw } );
$licence_text->windowCreate('text_center.last', -window => $button_OK_lic, -align => 'center');
$licence_window->protocol('WM_DELETE_WINDOW' => sub { $licence_window->withdraw } );
$licence_window->withdraw;
##Licence END

##Menu START
$mw->configure(-menu => my $menubar = $mw->Menu);
my $file_menu = $menubar->cascade(-label => '~File');
my $help_menu = $menubar->cascade(-label => '~Help');

my $new_working_directory = $file_menu->command(-label => 'New working directory', -underline => 0, -command => sub { new_working_directory(\$working_dir) } );
my $download_db = $file_menu->command(-label => 'Download enzyme database', -underline => 0, -command => \&download_enzyme_db);
my $save_selected_enzymes = $file_menu->command(-label => 'Save selected enzymes', -underline => 0, -command => sub { fileDialog_save_enzyme($mw) } );
my $load_selected_enzymes = $file_menu->command(-label => 'Load selected enzymes', -underline => 0, -command => sub { fileDialog_load_enzyme($mw) });
$file_menu->separator;
$file_menu->command(-label => 'Exit', -underline => 0, -command => sub { exit });

my $licence = $help_menu->command(-label => 'Licence information', -underline => 0, -command => sub { $licence_window->deiconify; $licence_window->raise });
my $about = $help_menu->command(-label => 'About vcf2caps', -underline => 0, -command => sub { $about_window->deiconify; $about_window->raise } );


##Menu END

my $folder_image = $mw->Photo(-file => 'icons/file_add.gif');
my $analyze_image = $mw->Photo(-file => 'icons/analyze.gif');
my $ok_image = $mw->Photo(-file => 'icons/ok.gif');
my $fail_image = $mw->Photo(-file => 'icons/fail.gif');
my $cancel = $mw->Photo(-file => 'icons/cancel.gif');

my $processing_gif = $mw->Animation(-format => 'gif', -file => 'icons/working.gif');

my $top = $mw->Frame->pack(-side => 'top', -fill => 'x');
my $bottom = $mw->Frame->pack(-side => 'top', -fill => 'both');

my $terminal = $bottom->Scrolled('Text', -scrollbars => 'e', -insertwidth => 0, -insertontime => 0, -insertofftime => 0, -wrap => 'word', -foreground => 'gray95' , -background => 'black', -relief => 'groove', -pady => 5, -padx => 5, -height => 15);
$terminal->tagConfigure('warning', -foreground => 'red');
$terminal->tagConfigure('warning_p', -foreground => 'red', -spacing3 => 0);
$terminal->tagConfigure('mark', -foreground => 'yellow');
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
#Scrolled
my $R_frame_enzymeProperties_listBox;

my $R_frame_allEnzymes_frame = $R_frame->Frame->pack(-side => 'left');
my $R_frame_buttons_frame = $R_frame->Frame->pack(-side => 'left');
my $R_frame_selEnzymes_frame = $R_frame->Frame->pack(-side => 'left');
$R_frame_allEnzymes_frame->Label(-text => 'All enzymes')->pack(-side => 'top',  -anchor => 'w', -padx => 10);
$R_frame_selEnzymes_frame->Label(-text => 'Selected enzymes')->pack(-side => 'top');
my $R_frame_selEnzymes_listBox;

my $R_frame_allEnzymes_listBox = $R_frame_allEnzymes_frame->Scrolled(
	'Listbox',
	-scrollbars => 'e',
	-height => 12,
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
#	-side => 'top',
	-anchor => 'c'
);
my $R_frame_selectAll_Button = $R_frame_buttons_frame->Button(
	-text => '>>',
	-width => 2,
	-command => sub {
		my $all_enz_index_size = $R_frame_allEnzymes_listBox->size();
#		my @all_enz_index = $R_frame_allEnzymes_listBox->selectionSet(0, 'end');
		for (my $i = 0; $i < $all_enz_index_size; $i++)
#		foreach my $index (@all_enz_index)
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
#	-side => 'top',
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
#	-side => 'top',
	-anchor => 'c'
);
$R_frame_selEnzymes_listBox = $R_frame_selEnzymes_frame->Scrolled(
	'Listbox',
	-scrollbars => 'e',
	-height => 12,
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
	#####################
	
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
	#			$R_frame_enzymeProperties_label_text->tagAdd('iso_1', '5.0', '5.13');
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
	#				$R_frame_enzymeProperties_label_text->tagAdd('iso_1', '5.0', '5.13');
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
		$R_frame_enzymeProperties_values_text->insert('end', join("", @seq) . "\n");
		$R_frame_enzymeProperties_values_text->insert('end', $data[0] . "\n");
		$R_frame_enzymeProperties_values_text->insert('end', $data[2]);
		
		$R_frame_enzymeProperties_label_text->tagConfigure('comp', -foreground => 'black');
	}
	#####################	
	
#	my @data = split("\t", $enzymes_db{$selected_enz_name});
#	my @seq;
#	for ( my $i = 0; $i < scalar( split("", $data[1]) ); $i++ )
#	{
#		if ($data[0] == $i) { push @seq, "'" }
#		else
#		{
#			push @seq, ( split("", $data[1]) )[$i]
#		}
#	}

#	$R_frame_enzymeProperties_values_text->delete('0.0', 'end');
#	$R_frame_enzymeProperties_values_text->insert('end', $selected_enz_name . "\n");
#	$R_frame_enzymeProperties_values_text->insert('end', join("", @seq) . "\n");
#	$R_frame_enzymeProperties_values_text->insert('end', $data[0] . "\n");
#	$R_frame_enzymeProperties_values_text->insert('end', $data[2]);
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
					if (!defined $x)
					{
						print join("\t", $enzymes_db{$enzyme_name});
						print "\n";
					}
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
####################################
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
			
#			$R_frame_enzymeProperties_label_text->tagAdd('iso_1', '5.0', '5.13');
		}
		else
		{
			my $first = 0;
			foreach my $isoschisomer ( split(",", ( split("\t", $enzymes_db{$selected_enz_name}) )[3] ) )
			{
				if ($first == 0) { $R_frame_enzymeProperties_values_text->insert('end', "$isoschisomer") }
				else { $R_frame_enzymeProperties_values_text->insert('end', "\n$isoschisomer") }
				$first++;
#				$R_frame_enzymeProperties_label_text->tagAdd('iso_1', '5.0', '5.13');
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
		$R_frame_enzymeProperties_values_text->insert('end', join("", @seq) . "\n");
		$R_frame_enzymeProperties_values_text->insert('end', $data[0] . "\n");
		$R_frame_enzymeProperties_values_text->insert('end', $data[2]);
		
		$R_frame_enzymeProperties_label_text->tagConfigure('iso', -foreground => 'black');
#		$R_frame_enzymeProperties_label_text->tagAdd('iso', '5.0', '5.13');
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
		$R_frame_enzymeProperties_values_text->insert('end', join("", @seq) . "\n");
		$R_frame_enzymeProperties_values_text->insert('end', $data[0] . "\n");
		$R_frame_enzymeProperties_values_text->insert('end', $data[2]);
		
		$R_frame_enzymeProperties_label_text->tagConfigure('comp', -foreground => 'black');
	}
});
$R_frame_enzymeProperties_values_text = $R_frame->Scrolled('Text', -scrollbars => 'e', -insertwidth => 0, -insertontime => 0, -insertofftime => 0, -width => 17, -wrap => 'word', -background => 'gray95', -relief => 'groove', -pady => 5, -padx => 5, -spacing3 => 10, -height => 6);
#$R_frame_enzymeProperties_values_text->insert('end', "Enzyme\n");
####################################

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

# $L_upper_1_frame - 1 oznacza nr wiersza
my $L_upper_1_frame = $L_upper_frame->Frame->pack(-side => 'top', -fill => 'x');
my $L_upper_2_frame = $L_upper_frame->Frame->pack(-side => 'top', -fill => 'x');
my $L_upper_3_frame = $L_upper_frame->Frame->pack(-side => 'top', -fill => 'x');
my $L_upper_4_frame = $L_upper_frame->Frame->pack(-side => 'top', -fill => 'x');
my $L_upper_5_frame = $L_upper_frame->Frame->pack(-side => 'top', -fill => 'x');

# $L_upper_1_2_frame - 1_2 - oznacza 2 kolumnę 1 wiersza
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
$L_upper_4_1_frame->Label(-text => 'Raw VCF', -width => 12)->pack(-side => 'top');
$L_upper_5_1_frame->Label(-text => 'Simplified VCF', -width => 12)->pack(-side => 'top');

$L_upper_1_2_frame->Label(-text => 'Selected file', -width => 23)->pack(-side => 'top');
my $enzyme_entry = $L_upper_2_2_frame->Entry(-insertwidth => 1, -width => 20,-textvariable => \$enzyme_file_name)->pack(-side => 'top',-pady => 1.3);
my $reference_entry = $L_upper_3_2_frame->Entry(-insertwidth => 1, -width => 20,-textvariable => \$reference_file_name)->pack(-side => 'top',-pady => 1.3);
my $raw_vcf_entry = $L_upper_4_2_frame->Entry(-insertwidth => 1, -width => 20,-textvariable => \$raw_vcf_file_name)->pack(-side => 'top',-pady => 1.3);
my $vcf_entry = $L_upper_5_2_frame->Entry(-insertwidth => 1, -width => 20,-textvariable => \$vcf_file_name)->pack(-side => 'top',-pady => 1.3);



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
my $L_lower_col2_mining_button;
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
		#	$enzyme_check_status->configure(-text => "Error - the file does not exist");
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
			#$reference_check_status->configure(-text => "Error - the file does not exist");
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
		if (defined $raw_vcf_file_name and -e $raw_vcf_file_name and $jobID == 0)
		{
			raw_start_vcf_check();
			$raw_vcf_analyze_button->configure(-state => 'disabled');
			$raw_chooseFile_button->configure(-state => 'disabled');
		}
		elsif (defined $raw_vcf_file_name and -e $raw_vcf_file_name and $jobID > 0)
		{
			$terminal->insert('end', "Warning", 'warning');
			$terminal->insert('end', " - the .\n");
			$terminal->see('end');
		}
		elsif (defined $raw_vcf_file_name and -e $raw_vcf_file_name and $reference_analysis_results[0] == 0)
		{
			$terminal->insert('end', "Warning", 'warning');
			$terminal->insert('end', " - before parsing VCF file, the reference file must first be loaded.\n\n");
			$terminal->see('end');
		}
		elsif (defined $raw_vcf_file_name and !-e $raw_vcf_file_name)
		{
			$vcf_check->configure(-image => $fail_image);
			#$vcf_check_status->configure(-text => "Error - the file does not exist");
			$terminal->insert('end', "Warning", 'warning');
			$terminal->insert('end', " - the file does not exist.\n\n");
			$terminal->see('end');
		}
		

	}
)->pack(
	-side => 'left'
);
my $vcf_chooseFile_button;
$vcf_chooseFile_button = $L_upper_5_3_frame->Button(
	-text => 'D',
	-image => $folder_image,
	-command => sub { fileDialog($mw,"vcf") }
)->pack(
	-side => 'left'
);
my $vcf_analyze_button;
$vcf_analyze_button = $L_upper_5_3_frame->Button(
	-text => 'D',
	-image => $analyze_image,
	-command => sub {
		if ($snps_seq_len !~ /^[0-9]+$/)
		{
			$terminal->insert('end', "Warning", 'warning');
			$terminal->insert('end', " - something is wrong with the marked parameter value. The parameter must be numerical.\n\n");
			$terminal->see('end');
			if ($snps_seq_len !~ /^[0-9]+$/) { $L_center_col2_1_entry->configure(-background => 'red') } else { $L_center_col2_1_entry->configure(-background => 'white') }
		}
		else
		{
			$L_center_col2_1_entry->configure(-background => 'white');
			
			
			if (defined $vcf_file_name and -e $vcf_file_name and $jobID == 0 and $reference_analysis_results[0] == 1)
			{
				start_vcf_check();
				$vcf_analyze_button->configure(-state => 'disabled');
				$vcf_chooseFile_button->configure(-state => 'disabled');
				$L_lower_col1_mining_button->configure(-state => 'disabled');
				$L_lower_col2_mining_button->configure(-state => 'disabled');
			}
			elsif (defined $vcf_file_name and -e $vcf_file_name and $jobID > 0 and $reference_analysis_results[0] == 1)
			{
				$terminal->insert('end', "Warning", 'warning');
				$terminal->insert('end', " - the .\n");
				$terminal->see('end');
			}
			elsif (defined $vcf_file_name and -e $vcf_file_name and $reference_analysis_results[0] == 0)
			{
				$terminal->insert('end', "Warning", 'warning');
				$terminal->insert('end', " - before parsing VCF file, the reference file must first be loaded.\n\n");
				$terminal->see('end');
			}
			elsif (defined $vcf_file_name and !-e $vcf_file_name)
			{
				$vcf_check->configure(-image => $fail_image);
				#$vcf_check_status->configure(-text => "Error - the file does not exist");
				$terminal->insert('end', "Warning", 'warning');
				$terminal->insert('end', " - the file does not exist.\n\n");
				$terminal->see('end');
			}
		}
		
		

	}
)->pack(
	-side => 'left'
);

$L_upper_1_4_frame->Label(
	-text => 'Status',
#	-width => 15
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

$vcf_check = $L_upper_5_4_frame->Label(
	-text => 'none',
	-foreground => 'grey',
)->pack(
	-side => 'left',
	-padx => 5
);
$vcf_check_status = $L_upper_5_4_frame->Label(
	-padx => 5
)->pack(
	-side => 'left'
);

#############

my $L_center_col1_frame = $L_center_frame->Frame->pack(-side => 'left', -fill => 'x', -pady => 10, -padx => 5);
my $L_center_col2_frame = $L_center_frame->Frame->pack(-side => 'left', -fill => 'x');
my $L_center_col1_1_label = $L_center_col1_frame->Label(-anchor => 'w', -text => "DNA sequence length flanking SNP/indel in the 'snps.txt'")->pack(-side => 'top', -fill => 'x');
my $L_center_col1_2_label = $L_center_col1_frame->Label(-anchor => 'w', -text => "DNA sequence length flanking SNP/indel in the output file")->pack(-side => 'top', -fill => 'x');
$L_center_col2_1_entry = $L_center_col2_frame->Entry(-insertwidth => 1, -width => 5, -textvariable => \$snps_seq_len, -justify => 'right')->pack(-side => 'top', -padx => 10);
$L_center_col2_2_entry = $L_center_col2_frame->Entry(-insertwidth => 1, -width => 5, -textvariable => \$output_seq_len, -justify => 'right')->pack(-side => 'top', -padx => 10);
my $L_lower_container_frame = $L_lower_frame->Frame->pack(-side => 'top', -fill => 'x');
my $L_lower_row0_frame = $L_lower_container_frame->Frame->pack(-side => 'top', -fill => 'x');
my $L_lower_row1_frame = $L_lower_container_frame->Frame->pack(-side => 'top', -fill => 'x');
my $L_lower_row2_frame = $L_lower_container_frame->Frame->pack(-side => 'top', -fill => 'x');

my $L_lower_col0_allSNPs = $L_lower_row0_frame->Checkbutton(
	-text => 'Mine CAPS from polymorphic SNPs/indels only',
	-variable => \$polymorphicSNPsOnly)->pack(-side => 'left', -anchor => 'w');

my $caps_singleCut_result_label;
$L_lower_col1_mining_button = $L_lower_row1_frame->Button(
	-text => 'Start CAPS mining',
	-width => 15,
	-state => 'disabled',
	-command => sub {
		if ($output_seq_len !~ /^[0-9]+$/ or ($custom == 1 and $custom_value !~ /^[0-9,]+$/) )
		{
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
				$caps_singleCut_result_label->packForget; # Usuwa etykietę z wynikami po filtrowaniu single-cut
				$L_lower_col1_mining_button->configure(-state => 'disabled');
				$fancy_SNPNo = 0;
				start_snp2caps();
			}
			elsif ( scalar(@selected_enz_names) == 0 )
			{
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
my $caps_mining_stop_button = $caps_mining_progress_frame->Button(-image => $cancel, -command => sub { $stop = 1 } )->pack(-side => 'left', -anchor => 'w');;
my $L_lower_col1_mining_textFrame = $caps_mining_progress_frame->Text(-width => 13, -height => 1, -state => 'disabled')->pack(-side => 'left', -anchor => 'w', -padx => 5);
my $progressBar = $L_lower_col1_mining_textFrame->ProgressBar(-variable => \$capsMining_percent, -width => 14, -length => 90, -gap => 0, -from => 0, -to => 100, -foreground => 'blue', -troughcolor => 'white');
$L_lower_col1_mining_textFrame->windowCreate('end', -window => $progressBar);
my $caps_mining_progress_label = $caps_mining_progress_frame->Label()->pack(-side => 'left', -anchor => 'w');


$L_lower_col2_mining_button = $L_lower_row2_frame->Button(
	-text => 'Select single-cut',
	-width => 15,
	-state => 'disabled',
	-command => sub {
		if ( scalar(@selected_enz_names) > 0 and $enzyme_analysis_results[0] == 1 and $reference_analysis_results[0] == 1 and $vcf_analysis_results{err_code} == 1 and $numberOfSNPs > 0 and $jobID == 0)
		{
			$L_lower_col2_mining_button->configure(-state => 'disabled');
			$numberOfSNPsBefore = 0;
			start_singleCutSite();
		}
	}
)->pack(-side => 'left', -anchor => 'w');
my $caps_singleCut_progress_frame = $L_lower_row2_frame->Frame;
$caps_singleCut_result_label = $L_lower_row2_frame->Label->pack(-side => 'left', -anchor => 'w');
my $caps_singleCut_stop_button = $caps_singleCut_progress_frame->Button(-image => $cancel, -command => sub { $stop = 1 } )->pack(-side => 'left', -anchor => 'w');;
my $caps_singleCut_textFrame = $caps_singleCut_progress_frame->Text(-width => 13, -height => 1, -state => 'disabled')->pack(-side => 'left', -anchor => 'w', -padx => 5);
my $singleCut_progressBar = $caps_singleCut_textFrame->ProgressBar(-variable => \$singleCutSite_percent, -width => 14, -length => 90, -gap => 0, -from => 0, -to => 100, -foreground => 'blue', -troughcolor => 'white');
$caps_singleCut_textFrame->windowCreate('end', -window => $singleCut_progressBar);
my $singleCut_progress_label = $caps_singleCut_progress_frame->Label()->pack(-side => 'left', -anchor => 'w');




$terminal->pack(-padx => 5, -pady => 5, -fill => 'x');
$terminal->insert('end', "vcf2CAPS v2.0\n\n");
$terminal->insert('end', "Welcome ...\n\n");


sub download_enzyme_db
{
	$terminal->insert('end', "Downloading the database from ");
	$terminal->insert('end', "http://rebase.neb.com/rebase/link_gcg", 'mark');
	$terminal->insert('end', " ...\n\n");
	$terminal->see('end');
	
	my $url = 'http://rebase.neb.com/rebase/link_gcg';
	my $file = $working_dir . "link_gcg";
	print "$file\n";
	my $response_code = getstore($url, $file);
	if ($response_code != 200)
	{
		$terminal->insert('end', "Warning", 'warning');
		$terminal->insert('end', " - something went wrong during download.\n\n");
		$terminal->see('end');
	}
	elsif ($response_code == 200)
	{
		$terminal->insert('end', "The database file downloaded sucessfully.\n\n");
		$terminal->see('end');
		
		$enzyme_file_name = $file;
		
		if (defined $enzyme_file_name and -f $file and $jobID == 0)
		{
			
			start_enzymes_check(); $enzyme_analyze_button->configure(-state => 'disabled');
		}
		elsif (defined $enzyme_file_name and !-f $file)
		{
			$enzyme_check->configure(-image => $fail_image);
		#	$enzyme_check_status->configure(-text => "Error - the file does not exist");
			$terminal->insert('end', "Warning", 'warning');
			$terminal->insert('end', " - the file does not exist.\n\n");
			$terminal->see('end');
		}
		
		
	}
}

sub new_working_directory
{
	my $working_dir_ref = $_[0];
	$$working_dir_ref = $mw->chooseDirectory(-initialdir => '.', -title => 'Choose a working directory');
	if (!defined $$working_dir_ref)
	{
		$terminal->insert('end', "No directory selected.\n\n");
		$terminal->see('end');
		$$working_dir_ref = "";
	}
	else
	{
		$terminal->insert('end', "Selected '");
		$terminal->insert('end', "$$working_dir_ref", 'mark');
		$terminal->insert('end', "' as a working directory.\n\n");
		$terminal->see('end');
		$$working_dir_ref = $$working_dir_ref . "/";
	}
	
}

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
				$terminal->insert('end', "The selected enzymes were saved sucessfully to the file '");
				$file =~ s/.*[\\\/]//g;
				$terminal->insert('end', "$file", 'mark');
				$terminal->insert('end', "'.\n\n");
				$terminal->see('end');
			}
			else
			{
				$terminal->insert('end', "Warning", 'warning');
				$terminal->insert('end', " - there was a problem during saving the file.\n\n");
				$terminal->see('end');
			}
		}
		else
		{
			$terminal->insert('end', "Warning", 'warning');
			$terminal->insert('end', " - the selected enzymes frame is empty. Please, select some enzymes and then try to save them.\n\n");
			$terminal->see('end');
		}
	}
}

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
			open my $fh, '<', $file or $err = 1;
				while (<$fh>)
				{
					chomp $_;
					if ($_ =~ /^[A-Za-z0-9-]+$/) { push @loaded_enzymes, $_ }
					else { $err = 2; last }
				}
			close $fh;
			
			if ($err == 1)
			{
				$terminal->insert('end', "Warning", 'warning');
				$terminal->insert('end', " - there was a problem during opening the file.\n\n");
				$terminal->see('end');
			}
			elsif ($err == 2)
			{
				$terminal->insert('end', "Warning", 'warning');
				$terminal->insert('end', " - something is wrong with the file content. Does it really contain the enzymes list?\n\n");
				$terminal->see('end');
			}

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
		$terminal->insert('end', "Warning", 'warning');
		$terminal->insert('end', " - please, load the enzyme database before loading selected enzymes.\n\n");
		$terminal->see('end');
	}
}

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

sub LOG
{
	my $text = $_[0];
	
	if ($log_first_use == 0)
	{
		if (-e "log.txt") { unlink "log.txt" }
		
		$log_first_use = 1;
		open my $Ofh, '>>', $working_dir . "log.txt";
			my $date = localtime();
			print $Ofh "vcf2caps v2.0\n";
			print $Ofh $date . "\n";
		close $Ofh;
	}
	
#	if ($^O eq 'MSWin32') { $new_line = "\r\n" }
#	elsif ($^O eq 'darwin') { $new_line = "\r" }
#	else { $new_line = "\n" }
	
	open my $Ofh, '>>', $working_dir . "log.txt";
		print $Ofh $text . "\n";
	close $Ofh;
}

MainLoop;



sub start_enzymes_check
{
	if (defined $enzyme_file_name and -e $enzyme_file_name)
	{
		@enzyme_analysis_results = (0);
		$jobID = 2;
	#	reference_check();
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
				#	$enzyme_check_status->configure(-text => "OK");
#					$L_upper_2_4_frame->Label(-text => "OK")->pack(-side => 'left');
					
					$R_frame_allEnzymes_listBox->delete(0,'end');
					
					foreach my $enzyme_name (@allEnzymesNames)
					{
						$R_frame_allEnzymes_listBox->insert('end', $enzyme_name);
					}
					
					$R_frame_allEnzymes_No->configure(-text => scalar(@allEnzymesNames) . " enzymes" );
					$terminal->insert('end', "Loaded REBASE database (version $enzymes_db{date}) comprising " . scalar(@allEnzymesNames) . " enzymes.\n\n");
					$terminal->see('end');
					
					if ($reference_analysis_results[0] == 1 and $enzyme_analysis_results[0] == 1 and $vcf_analysis_results{err_code} == 1)
					{
						$L_lower_col1_mining_button->configure(-state => 'normal');
					}
				}
				elsif ($enzyme_analysis_results[0] == 3)
				{
					$enzyme_check->configure(-image => $fail_image);
				#	$enzyme_check_status->configure(-text => "Error - something is wrong with the database file.");
					$terminal->insert('end', "Warning", 'warning');
					$terminal->insert('end', " - something is wrong with the database file.\n\n");
					$terminal->see('end');
				}
				
				
#				print "$enzymes_db{date}\n";
#				foreach my $company ( split("\t", $enzymes_db{companies}) )
#				{
#					my @data = split(",", $company);
#					print "$data[0] -> $data[1]\n";
#				}
#				print "EcoRI:\n";
				

#					my @data = split(",", $enzymes_db{EcoRI});
#					print "Cut site: $data[0]\n";
#					print "Sequence: $data[1]\n";
#					print "Company: $data[4]\n";
					

				$enzyme_analyze_button->configure(-state => 'normal');
				$enzyme_chooseFile_button->configure(-state => 'normal');
				$repeat->cancel;
			}
	#		print "loop: $c\n";
		} );
		
	}
	elsif (defined $enzyme_file_name)
	{
		$enzyme_check->configure(-image => $fail_image);
		#$enzyme_check_status->configure(-text => "Error - the file does not exist");
		$terminal->insert('end', "Warning", 'warning');
		$terminal->insert('end', " - the file does not exist.\n\n");
		$terminal->see('end');
	}
	
}

sub start_reference_check
{
	if (defined $reference_file_name and -e $reference_file_name)
	{
		@reference_analysis_results = (0);
		$jobID = 1;
		my $reference_file_name_tmp;
		$reference_file_name_tmp = $reference_file_name;
		$reference_file_name_tmp =~ s/\..+$//g;
		if (-e $reference_file_name_tmp . ".index") { unlink $reference_file_name_tmp . ".index" }
		$reference_file_name_tmp = $reference_file_name;
		$reference_file_name_tmp =~ s/.*[\\\/]//g;
		
		$terminal->insert('end', "Start reference file '");
		$terminal->insert('end', "$reference_file_name_tmp", 'mark');
		$terminal->insert('end', "' integrity check ...\n\n");
		$terminal->see('end');
		
	#	reference_check();
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
				#	$reference_check_status->configure(-text => "OK");
				
					
					$terminal->insert('end', "Integrity of the reference file '");
					$terminal->insert('end', "$reference_file_name_tmp", 'mark');
					$terminal->insert('end', "' confirmed. The index file '");
					$reference_file_name_tmp =~ s/(?<=\.).+$//;
					$terminal->insert('end', "$reference_file_name_tmp" . "index", 'mark');
					$terminal->insert('end', "' was created.\n\n");
					$terminal->see('end');
					
					if ($samtools_error_code[1] == 1 and $reference_analysis_results[0] == 1 and $enzyme_analysis_results[0] == 1 and $vcf_analysis_results{err_code} == 1)
					{
						$L_lower_col1_mining_button->configure(-state => 'normal');
					}
				}
				elsif ($reference_analysis_results[0] == 3)
				{
					$reference_check->configure(-image => $fail_image);
				#	$reference_check_status->configure(-text => "Error - duplicated sequence: $reference_analysis_results[1]");
					$terminal->insert('end', "Warning", 'warning');
					$terminal->insert('end', " - duplicated sequence: $reference_analysis_results[1].\n\n");
					$terminal->see('end');
				}
				elsif ($reference_analysis_results[0] == 4)
				{
					$reference_check->configure(-image => $fail_image);
					my @err_data = split(",", $reference_analysis_results[1]);
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
				
				$reference_analyze_button->configure(-state => 'normal');
				$reference_chooseFile_button->configure(-state => 'normal');
				$repeat->cancel;
			}
	#		print "loop: $c\n";
		} );
	}
	elsif (defined $reference_file_name)
	{
		$reference_check->configure(-image => $fail_image);
		#$reference_check_status->configure(-text => "Error - the file does not exist");
		
		$terminal->insert('end', "Warning", 'warning');
		$terminal->insert('end', " - the file does not exist.\n\n");
		$terminal->see('end');
	}
}

sub raw_start_vcf_check
{
	$raw_vcf_analysis_results{err_code} = 0;
	$jobID = 7;
	my $repeat;
	$raw_vcf_check->configure(-image => $processing_gif);
	$processing_gif->start_animation;
	
	my $raw_vcf_file_name_tmp = $raw_vcf_file_name;
	$raw_vcf_file_name_tmp =~ s/.*[\\\/]//g;
	my $raw_vcf_tmp = $raw_vcf_file_name_tmp;
	$raw_vcf_tmp =~ s/\..+$//g;
	$raw_vcf_tmp = $raw_vcf_tmp . ".svcf";
	
	if (-e $raw_vcf_tmp) { unlink $raw_vcf_tmp }
	
	$terminal->insert('end', "Start raw VCF file '");
	$terminal->insert('end', "$raw_vcf_file_name_tmp", 'mark');
	$terminal->insert('end', "' convertion to the simplified form ...\n\n");
	$terminal->see('end');
	
	$repeat = $mw->repeat( 100 => sub {
		if ($raw_vcf_analysis_results{err_code} != 0)
		{
			$processing_gif->stop_animation;
			
			if ($raw_vcf_analysis_results{err_code} == 1)
			{
				$raw_vcf_check->configure(-image => $ok_image);
				
				$terminal->insert('end', "Convertion of the VCF file '");
				$terminal->insert('end', "$raw_vcf_file_name_tmp", 'mark');
				$terminal->insert('end', "' completed. The simplified form was saved in: '");
				
				$terminal->insert('end', "$working_dir" . "$raw_vcf_tmp", 'mark');
				$terminal->insert('end', "'.\n\n");
				$terminal->see('end');
				
				
			}
			elsif ($raw_vcf_analysis_results{err_code} == 2)
			{
				$raw_vcf_check->configure(-image => $fail_image);

				$terminal->insert('end', "Warning", 'warning');
				$terminal->insert('end', " - the file '");
				$terminal->insert('end', "$raw_vcf_file_name_tmp", 'mark');
				$terminal->insert('end', "' does not have the header '");
				$terminal->insert('end', "#CHROM",'mark');
				$terminal->insert('end', "'. Is it really VCF file?\n\n");
				$terminal->see('end');
			}
			elsif ($raw_vcf_analysis_results{err_code} == 3)
			{
				$raw_vcf_check->configure(-image => $fail_image);

				$terminal->insert('end', "Warning", 'warning');
				$terminal->insert('end', " - something is wrong with the file '");
				$terminal->insert('end', "$raw_vcf_file_name_tmp", 'mark');
				$terminal->insert('end', "'.\n\n");
				$terminal->see('end');
			}
			
			$raw_vcf_analyze_button->configure(-state => 'normal');
			$raw_chooseFile_button->configure(-state => 'normal');
			$repeat->cancel;
		}
	});
}

sub start_vcf_check
{
	%vcf_analysis_results = ();
	%vcf_analysis_results = (err_code => 0);
	$jobID = 3;
	my $repeat;
	my $c = 0;
	$vcf_check->configure(-image => $processing_gif);
	$processing_gif->start_animation;
	my $vcf_file_name_tmp = $vcf_file_name;
	$vcf_file_name_tmp =~ s/.*[\\\/]//g;
	$terminal->insert('end', "Start VCF file '");
	$terminal->insert('end', "$vcf_file_name_tmp", 'mark');
	$terminal->insert('end', "' integrity check ...\n\n");
	$terminal->see('end');
	my $next_step = 0;
	$repeat = $mw->repeat( 100 => sub {
		if ($vcf_analysis_results{err_code} != 0)
		{
			#$processing_gif->stop_animation;
			
			if ($vcf_analysis_results{err_code} == 1)
			{
			#	$vcf_check->configure(-image => $ok_image);
			#	$vcf_check_status->configure(-text => "OK - No. of individuals: $vcf_analysis_results{NoOfIndv}, No. of SNPs: $vcf_analysis_results{NoOfSNPs}");
				
				$terminal->insert('end', "Integrity check of the VCF file '");
				$terminal->insert('end', "$vcf_file_name_tmp", 'mark',);
				$terminal->insert('end', "' confirmed. Statistics:\n- number of individuals: ");
				$terminal->insert('end', "$vcf_analysis_results{NoOfIndv}", 'mark');
				$terminal->insert('end', "\n- number of SNPs/indels: ",);
				$terminal->insert('end', "$vcf_analysis_results{NoOfSNPs}\n\n", 'mark');
				$terminal->see('end');
				
				$terminal->insert('end', "Parsing the '");
				$terminal->insert('end', "$vcf_file_name_tmp", 'mark');
				$terminal->insert('end', "' file ... ");
				$terminal->insert('end', "0.0%", 'percent');
				$terminal->see('end');
				$next_step = 1;
				$jobID = 4;
				
				
				my $percent_index_start;
				$percent_index_start = $terminal->search(-regexp, -backwards => '[0-9]+\.[0-9]%', 'end');
				my $repeat2;
				$repeat2 = $mw->repeat( 100 => sub {
					if ($samtools_error_code[1] != 0)
					{
						print "\$samtools_error_code[1]: $samtools_error_code[1]\n";
					#	if (defined $percent_index_start)
					#	{
							my $percent = ( ($line_vcf - 1) / $vcf_analysis_results{NoOfSNPs}) * 100;
							my $percent_index_stop = $terminal->search(-regexp, -backwards => '%', 'end');
							$percent_index_stop += 1;
							$terminal->delete("$percent_index_start", "$percent_index_stop");
							$terminal->insert("$percent_index_start", sprintf ("%.1f%%", $percent));
					#	}
						
						$processing_gif->stop_animation;
						
						my $vcf_file_name_tmp = $vcf_file_name;
						$vcf_file_name_tmp =~ s/.*[\\\/]//g;
						my $reference_file_name_tmp = $reference_file_name;
						$reference_file_name_tmp =~ s/.*[\\\/]//g;
						
						my $sequencesNotPresentInRef_No = scalar(@sequencesNotPresentInRef);
						my $sequencesNotPresentInRef_No_forReport = scalar(@sequencesNotPresentInRef);
						my $markersOnTheEdge_No = scalar(@markersOnTheEdge);
						my $markersOnTheEdge_No_forReport = scalar(@markersOnTheEdge);
						
						if ($samtools_error_code[1] == 1 and $sequencesNotPresentInRef_No == 0 and $markersOnTheEdge_No == 0 )
						{
							$vcf_check->configure(-image => $ok_image);
							
							$terminal->insert('end', "\n\nParsing finished sucessfully.\n\n");
							$terminal->see('end');
							
							if ($samtools_error_code[1] == 1 and $reference_analysis_results[0] == 1 and $enzyme_analysis_results[0] == 1 and $vcf_analysis_results{err_code} == 1)
							{
								$L_lower_col1_mining_button->configure(-state => 'normal');
							}
						}
						elsif ($samtools_error_code[1] == 1 and $snpsNo > 0 and ($sequencesNotPresentInRef_No > 0 or $markersOnTheEdge_No > 0) )
						{
							
						
							$vcf_check->configure(-image => $ok_image);
							
							
							$terminal->insert('end', "\n\n");
							$terminal->insert('end', "Warning",'warning');
							if ( $sequencesNotPresentInRef_No > 0 and $markersOnTheEdge_No > 0)
							{
								$terminal->insert('end', " - parsing finished. However, some problems occurred:\n");
							}
							else
							{
								$terminal->insert('end', " - parsing finished. However, a problem occurred:\n");
							}
							
							if ( $sequencesNotPresentInRef_No > 0 )
							{
								if ( $sequencesNotPresentInRef_No <= 5 )
								{
									$terminal->insert('end', "Below listed SNPs/indels were located in the sequences not present in the reference file '");
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
									LOG("\n# VCF to sVCF conversion");
									LOG("# SNPs/indels located in the sequences not present in the reference file $reference_file_name_tmp");
									foreach my $text (@sequencesNotPresentInRef)
									{
										LOG($text);
									}
									$terminal->insert('end', "$sequencesNotPresentInRef_No_forReport", 'mark');
									$terminal->insert('end', " SNPs/indels were located in the sequences not present in the reference file '");
									$terminal->insert('end', "$reference_file_name_tmp", 'mark');
									$terminal->insert('end', "'. The list of those SNPs/indels were saved in the ");
									$terminal->insert('end', "log.txt", 'mark');
									$terminal->insert('end', " file.\n");
									$terminal->see('end');
								}
							}
							
							if ( $markersOnTheEdge_No > 0)
							{															
								if ( $markersOnTheEdge_No <= 5 )
								{
									$terminal->insert('end', "Below listed SNPs/indels were located closer to the edges of the sequences than ");
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
									LOG("\n# VCF to sVCF conversion");
									LOG("# SNPs/indels located closer to the edges of the sequences than 40 bp");
									foreach my $text (@markersOnTheEdge)
									{
										LOG($text);
									}
									$terminal->insert('end', "$markersOnTheEdge_No_forReport", 'mark');
									$terminal->insert('end', " SNPs/indels were located closer to the edges of the sequences than ");
									$terminal->insert('end', "$snps_seq_len", 'mark');
									$terminal->insert('end', " bp:\n");
									$terminal->insert('end', "The list of those SNPs/indels were saved in the ");
									$terminal->insert('end', "log.txt", 'mark');
									$terminal->insert('end', " file.\n");
									$terminal->see('end');
								}
							}
							
							$terminal->insert('end', "\n");
							$terminal->see('end');
							
							if ($samtools_error_code[1] == 1 and $reference_analysis_results[0] == 1 and $enzyme_analysis_results[0] == 1 and $vcf_analysis_results{err_code} == 1)
							{
								$L_lower_col1_mining_button->configure(-state => 'normal');
							}	
						}
						elsif ($snpsNo == 0)
						{
							$vcf_check->configure(-image => $fail_image);
							
							$terminal->insert('end', "\n\n");
							$terminal->insert('end', "Warning",'warning');
							$terminal->insert('end', " - parsing failed. Something went horribly wrong. Please, check the VCF file '");
							$terminal->insert('end', "$vcf_file_name_tmp",'mark');
							$terminal->insert('end', "' for any issues with the file format.\n\n");
							$terminal->see('end');
						}
						
				#		$vcf_analyze_button->configure(-state => 'normal');
				#		$vcf_chooseFile_button->configure(-state => 'normal');
						
						$samtools_error_code[1] = 0;
						$line_vcf = 0;
						$repeat2->cancel;
						
					}
					#print "\$line_vcf: $line_vcf\t\$vcf_analysis_results{NoOfSNPs}: $vcf_analysis_results{NoOfSNPs}\n";
					elsif ($line_vcf > 0 and $vcf_analysis_results{NoOfSNPs} > 0)
					{
						my $percent = ( ($line_vcf - 1) / $vcf_analysis_results{NoOfSNPs}) * 100;
						
						
						if (!defined $percent_index_start)
						{
							$percent_index_start = $terminal->search(-regexp, -backwards => '[0-9]+\.[0-9]%', 'end');
						}
						#print "$percent_index_start\n";
						my $percent_index_stop = $terminal->search(-regexp, -backwards => '%', 'end');
						$percent_index_stop += 1;
						$terminal->delete("$percent_index_start", "$percent_index_stop");
						#$terminal->see("$percent_index");
						$terminal->insert("$percent_index_start", sprintf ("%.1f%%", $percent));
					}
				});
			}
			elsif ($vcf_analysis_results{err_code} == 3)
			{
				$vcf_check->configure(-image => $fail_image);
			#	$vcf_check_status->configure(-text => "Error - the file does not have a header");
				
				
				$terminal->insert('end', "Warning", 'warning');
				$terminal->insert('end', " - the file '");
				$terminal->insert('end', "$vcf_file_name_tmp", 'mark');
				$terminal->insert('end', "' does not have the header '");
				$terminal->insert('end', "#CHROM",'mark');
				$terminal->insert('end', "'. Is it really sVCF file?\n\n");
				$terminal->see('end');
			}
			

			$vcf_analyze_button->configure(-state => 'normal');
			$vcf_chooseFile_button->configure(-state => 'normal');
			$repeat->cancel;
		}

	} );
}

sub start_snp2caps
{
	$jobID = 5;
	$L_lower_col2_mining_button->configure(-state => 'disabled');
	$capsMining_percent = 0;
	$enzyme_zip = 0;
	$enzyme_zip_stop = 0;
#	$progressBar->configure(-value => 0);
	$caps_mining_result_label->packForget;
	my $repeat;
	my $enzymes_for_analysis_No = keys %enzymes_for_analysis;
	
	$repeat = $mw->repeat( 100 => sub {
		$capsMining_percent = ( ($enzyme_zip) / $enzymes_for_analysis_No) * 100;
		if ($enzyme_zip > 0 and $capsMining_percent < 100)
		{
			$caps_mining_prepare_enzymes_label->pack(-side => 'left', -anchor => 'w', -padx => 5);
			
			$caps_mining_prepare_enzymes_label->configure(-text => sprintf ("Preparing enzymes list   %.1f%%", $capsMining_percent) );
		}
		elsif ($capsMining_percent == 100)
		{
			$caps_mining_prepare_enzymes_label->pack(-side => 'left', -anchor => 'w', -padx => 5);
			$caps_mining_prepare_enzymes_label->configure(-text => sprintf ("Preparing enzymes list   %.1f%%", $capsMining_percent) );
			$caps_mining_prepare_enzymes_label->packForget;
			$capsMining_percent = 0;
			$caps_mining_progress_label->configure(-text => sprintf ("%d/%d   %.1f%%", $fancy_SNPNo,$vcf_analysis_results{NoOfSNPs},$capsMining_percent) );
			$caps_mining_progress_frame->pack(-side => 'left', -anchor => 'w', -padx => 5);
			
			$snp2caps_results[1] = 0;
			$terminal->insert('end', "Start CAPS mining ...\n\n");
			$terminal->see('end');
			my $repeat2;
			$repeat2 = $mw->repeat( 100 => sub {
				
				
				
				if ($snp2caps_results[1] != 0)
				{
					
					
					@markersOnTheEdge = uniq(@markersOnTheEdge);
					
					my $markersOnTheEdge_No = scalar(@markersOnTheEdge);
					my $markersOnTheEdge_No_forReport = scalar(@markersOnTheEdge);
					
					if ($snp2caps_results[1] == 1 and $markersOnTheEdge_No == 0)
					{
						$caps_mining_progress_frame->packForget;
						my $numberOfSNPs_tmp = $numberOfSNPs;
						$numberOfSNPs_tmp =~ s/(\d{1,3}?)(?=(\d{3})+$)/$1 /g;
						$caps_mining_result_label->configure(-text => "$numberOfSNPs_tmp CAPS found");
						$caps_mining_result_label->pack(-side => 'left', -anchor => 'w');
												
						$terminal->insert('end', "CAPS mining finished.");
						$terminal->insert('end', " $numberOfSNPs_tmp", 'mark');
						$terminal->insert('end', " CAPS were saved to file '");
						$terminal->insert('end', "out.txt", 'mark');
						$terminal->insert('end', "'\n\n");;
						$terminal->see('end');
						
						$L_lower_col2_mining_button->configure(-state => 'normal');
					}
					elsif ($snp2caps_results[1] == 1 and $markersOnTheEdge_No > 0)
					{
						$caps_mining_progress_frame->packForget;
						my $numberOfSNPs_tmp = $numberOfSNPs;
						$numberOfSNPs_tmp =~ s/(\d{1,3}?)(?=(\d{3})+$)/$1 /g;
						$caps_mining_result_label->configure(-text => "$numberOfSNPs_tmp CAPS found");
						$caps_mining_result_label->pack(-side => 'left', -anchor => 'w');
						
						my $vcf_file_name_tmp = $vcf_file_name;
						$vcf_file_name_tmp =~ s/.*[\\\/]//g;
						
						$terminal->insert('end', "Warning",'warning');
						$terminal->insert('end', " - CAPS mining finished.");
						$terminal->insert('end', " $numberOfSNPs_tmp", 'mark');
						$terminal->insert('end', " CAPS were saved to file '");
						$terminal->insert('end', "out.txt", 'mark');
						$terminal->insert('end', "'. However, a problem occurred:\n");
												
						if ( $markersOnTheEdge_No <= 5 )
						{
							$terminal->insert('end', "Below listed SNPs/indels were located closer to the edges of the sequences than ");
							$terminal->insert('end', "$output_seq_len", 'mark');
							$terminal->insert('end', " bp:\n");
						
							foreach my $markerOnTheEdge (@markersOnTheEdge)
							{
								$terminal->insert('end', "- ");
								$terminal->insert('end', "$markerOnTheEdge\n\n",'mark');
							}
							$terminal->see('end');
						}
						else
						{
							LOG("\n# CAPS mining");
							LOG("# SNPs/indels located closer to the edges of the sequences than 500 bp");
							foreach my $text (@markersOnTheEdge)
							{
								LOG($text);
							}
							$terminal->insert('end', "$markersOnTheEdge_No_forReport", 'mark');
							$terminal->insert('end', " SNPs/indels were located closer to the edges of the sequences than ");
							$terminal->insert('end', "$output_seq_len", 'mark');
							$terminal->insert('end', " bp:\n");
							$terminal->insert('end', "The list of those SNPs/indels were saved in the ");
							$terminal->insert('end', "log.txt", 'mark');
							$terminal->insert('end', " file.\n");
							$terminal->see('end');
						}

						$terminal->insert('end', "\n");
						$terminal->see('end');
						
						$L_lower_col2_mining_button->configure(-state => 'normal');
					}
					elsif ($snp2caps_results[1] == 2)
					{
						$caps_mining_progress_frame->packForget;
						my $numberOfSNPs_tmp = $numberOfSNPs;
						$numberOfSNPs_tmp =~ s/(\d{1,3}?)(?=(\d{3})+$)/$1 /g;
						$caps_mining_result_label->configure(-text => "$numberOfSNPs_tmp CAPS found");
						$caps_mining_result_label->pack(-side => 'left', -anchor => 'w');
						
						$terminal->insert('end', "Warning", 'warning');
						$terminal->insert('end', " - CAPS mining canceled.");
						$terminal->insert('end', " $numberOfSNPs_tmp", 'mark');
						$terminal->insert('end', " CAPS were saved to file '");
						$terminal->insert('end', "out.txt", 'mark');
						$terminal->insert('end', "'\n\n");;
						$terminal->see('end');
						
						$L_lower_col2_mining_button->configure(-state => 'normal');
					}
					elsif ($snp2caps_results[1] == 2 and $markersOnTheEdge_No > 0)
					{
						$caps_mining_progress_frame->packForget;
						my $numberOfSNPs_tmp = $numberOfSNPs;
						$numberOfSNPs_tmp =~ s/(\d{1,3}?)(?=(\d{3})+$)/$1 /g;
						$caps_mining_result_label->configure(-text => "$numberOfSNPs_tmp CAPS found");
						$caps_mining_result_label->pack(-side => 'left', -anchor => 'w');
						
						my $vcf_file_name_tmp = $vcf_file_name;
						$vcf_file_name_tmp =~ s/.*[\\\/]//g;
						$terminal->insert('end', "Warning", 'warning');
						$terminal->insert('end', " - CAPS mining canceled.");
						$terminal->insert('end', " $numberOfSNPs_tmp", 'mark');
						$terminal->insert('end', " CAPS were saved to file '");
						$terminal->insert('end', "out.txt", 'mark');
						$terminal->insert('end', "'. However, a problem occurred:\n");
						
						if ( $markersOnTheEdge_No <= 5 )
						{
							$terminal->insert('end', "Below listed SNPs/indels were located closer to the edges of the sequences than ");
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
							LOG("# SNPs/indels located closer to the edges of the sequences than 500 bp");
							foreach my $text (@markersOnTheEdge)
							{
								LOG($text);
							}
							$terminal->insert('end', "$markersOnTheEdge_No_forReport", 'mark');
							$terminal->insert('end', " SNPs/indels were located closer to the edges of the sequences than ");
							$terminal->insert('end', "$snps_seq_len", 'mark');
							$terminal->insert('end', " bp:\n");
							$terminal->insert('end', "The list of those SNPs/indels were saved in the ");
							$terminal->insert('end', "log.txt", 'mark');
							$terminal->insert('end', " file.\n");
							$terminal->see('end');
						}
						
						
						$terminal->insert('end', "\n");
						$terminal->see('end');
					}
					
					
					
					$L_lower_col1_mining_button->configure(-state => 'normal');
					
					
					$repeat2->cancel;
				}
				elsif ($fancy_SNPNo > 0 and $vcf_analysis_results{NoOfSNPs} > 0)
				{
					$capsMining_percent = ( ($fancy_SNPNo) / $vcf_analysis_results{NoOfSNPs}) * 100;
					
					$caps_mining_progress_label->configure(-text => sprintf ("%d/%d   %.1f%%", $fancy_SNPNo,$vcf_analysis_results{NoOfSNPs},$capsMining_percent) );

					
					#$terminal->insert("$percent_index_start", sprintf ("%.1f%%", $percent));
				}
			});
			
			$repeat->cancel;
		}
		
	});
	
	
	
	
	
}

sub start_singleCutSite
{
	$jobID = 6;
	$singleCutSite_percent = 0;
	$caps_singleCut_result_label->packForget;
	$singleCut_progress_label->configure(-text => sprintf ("%d/%d   %.1f%%", $numberOfSNPsBefore,$numberOfSNPs,$singleCutSite_percent) );
	$caps_singleCut_progress_frame->pack(-side => 'left', -anchor => 'w', -padx => 5);
	my $repeat;
	$singleCutSite_results[1] = 0;
	$terminal->insert('end', "Start single-cut CAPS filtering ...\n\n");
	$terminal->see('end');
	$repeat = $mw->repeat( 100 => sub {
		if ($singleCutSite_results[1] != 0)
		{
			if ($singleCutSite_results[1] == 1)
			{
				$caps_singleCut_progress_frame->packForget;
				$numberOfSNPsAfter =~ s/(\d{1,3}?)(?=(\d{3})+$)/$1 /g;
				$caps_singleCut_result_label->configure(-text => "$numberOfSNPsAfter CAPS found with single-cut site");
				$caps_singleCut_result_label->pack(-side => 'left', -anchor => 'w');
				$terminal->insert('end', "Single-cut CAPS filtering finished.");
				$terminal->insert('end', " $numberOfSNPsAfter", 'mark');
				$terminal->insert('end', " CAPS were saved to '");
				$terminal->insert('end', "out_single.txt", 'mark');
				$terminal->insert('end', "' file\n\n");
				$terminal->see('end');
			}
			elsif ($singleCutSite_results[1] == 2)
			{
				$caps_singleCut_progress_frame->packForget;
				$numberOfSNPsAfter =~ s/(\d{1,3}?)(?=(\d{3})+$)/$1 /g;
				$caps_singleCut_result_label->configure(-text => "$numberOfSNPsAfter CAPS found with single-cut site");
				$caps_singleCut_result_label->pack(-side => 'left', -anchor => 'w');
				$terminal->insert('end', "Warning",'warning');
				$terminal->insert('end', " - single-cut CAPS filtering canceled.");
				$terminal->insert('end', " $numberOfSNPsAfter", 'mark');
				$terminal->insert('end', " CAPS were saved to '");
				$terminal->insert('end', "out_single.txt", 'mark');
				$terminal->insert('end', "' file\n\n");
				$terminal->see('end');
			}
			$L_lower_col2_mining_button->configure(-state => 'normal');
			
			
			$repeat->cancel;
		}
		
		elsif ($numberOfSNPsBefore > 0 and $vcf_analysis_results{NoOfSNPs} > 0)
		{
			$singleCutSite_percent = ( ($numberOfSNPsBefore) / $numberOfSNPs ) * 100;

			$singleCut_progress_label->configure(-text => sprintf ("%d/%d   %.1f%%", $numberOfSNPsBefore,$numberOfSNPs,$singleCutSite_percent) );
			
			#$terminal->insert("$percent_index_start", sprintf ("%.1f%%", $percent));
		}
	});
}

sub fileDialog {
	my $w = shift;
	my $file_type = shift;
	#    my $ent = shift;
	my $types;
	my $file;
	my @types_reference =
		(["Fasta files", [qw/.fa .fasta/],'TEXT'],
		["Text files", '.txt','TEXT'],
		["All files", '*']
	);
	my @types_vcf =
		(["VCF files", '.vcf', 'TEXT'],
		["Text files", '.txt'],
		["All files", '*']
	);
	my @types_svcf =
		(["Simplified VCF files", '.svcf', 'TEXT'],
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
			
			$reference_entry->delete(0, 'end');
			$reference_entry->insert(0, $file);
			$reference_entry->xview('end');
		}
	}
	elsif ($file_type eq "vcf")
	{
		$file = $w->getOpenFile(-filetypes => \@types_svcf, -defaultextension => '.svcf');
		if (defined $file and $file ne "")
		{
			$vcf_check->configure(-image => '', -text => 'none');
			
			$vcf_entry->delete(0, 'end');
			$vcf_entry->insert(0, $file);
			$vcf_entry->xview('end');
		}
	}
	elsif ($file_type eq "raw_vcf")
	{
		$file = $w->getOpenFile(-filetypes => \@types_vcf, -defaultextension => '.vcf');
		if (defined $file and $file ne "")
		{
			$raw_vcf_check->configure(-image => '', -text => 'none');
			
			$raw_vcf_entry->delete(0, 'end');
			$raw_vcf_entry->insert(0, $file);
			$raw_vcf_entry->xview('end');
		}
	}

	#    if (defined $file and $file ne '') {
	#	$ent->delete(0, 'end');
	#	$ent->insert(0, $file);
	#	$ent->xview('end');
	#   }
}



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
					my $prefix = $reference_file_name;
					$prefix =~ s/\.fa[a-z]*$//;
					my $line_len = 0; # Długość sekwencji nukleotydowej w każdej linii (za wyjątkiem ostatniej linii danej sekwencji) pliku fasta
					my $chrom_ID = ""; # Nazwa sekwencji z pliku referencji
					my $chrom_len = 0; # Całkowita długość danej sekwencji
					my $last_line_len = 0; # Długość ostatniej linii danej sekwencji
					my $alert = 0; # Jeżeli 1 - różna długość linii z sekwencjami w pliku fasta. Jeżeli tylko ostatnia linia chromosomu ma inną długość to ten parametr jest zerowany
					my @all_chrom_ID = (); # Zmienna przechowuje nazwy wszystkich sekwencji obecnych w referencji
					open my $fh, "<", $reference_file_name or die "Cannot open the file '$reference_file_name'.";
					open my $Ofh, '>>', $prefix . ".index" or die "Cannot create the file '$prefix.index'\n";
					
						L: while (<$fh>)
						{
					###########################
						
							chomp $_;
							if ($_ =~ /^>/) # Jeżeli linia rozpoczyna się od znaku zachęty '>', to ...
							{
								if ($chrom_len != 0) {print $Ofh "$last_line_len\t$chrom_len\n"; $chrom_len = 0} # ... sprawdza, czy zmienna $chrom_len jest różna od zera. Jeżeli tak, to do pliku indeksu dopisywana jest: długość ostatniej linii danej sekwencji oraz długość całej sekwencji. Następnie zmienna $chrom_len jest zerowana. Ta linijka kodu jest wykonywana, aby zapisać do pliku indeksu zebrane wcześniej (z poprzdenij rundy pętli) dane dla poprzedniej sekencji
								$chrom_ID = $_; # Zmiennej $chrom_ID przypisywana jest nazwa danej sekwencji
								if (scalar(@all_chrom_ID) > 0) # Jeżeli tablica @all_chrom_ID jest niezerowa (są zapisane w niej nazwy jakichś sekwencji), to ...
								{
									foreach my $c (@all_chrom_ID) # ... dla każdy element z tej tablicy ...
									{
										if ($chrom_ID eq $c) # ... sprawdza, czy jest on równy nazwie aktualnie badanej sekwencji. Jeżeli tak (jest duplikacja nazw sekwencji w pliku referencji), to ...
										{
				#									print "[\e[91mfail\e[0m]\n\n";
				#									print "Sequence name '$chrom_ID' duplicate present.\nExiting ...\n";
											unlink $prefix . ".index"; # ... usuwa utworzony plik indeksu, następnie ...
											
											@genome_error = (1,$chrom_ID);
											last L;
										}
									}
								}
								push @all_chrom_ID, $chrom_ID; # Po weryfikacji duplikacji nazw sekwencji do tablicy @all_chrom_ID dodawana jest nazwa aktualnie analizowanej sekwencji, a następnie ...
								print $Ofh $_ . "\t" . tell($fh) . "\t"; # ... do pliku indeksu jest zapisywana nazwa aktualnie analizowanej sekwencji oraz odpowiadająca jej pozycja kursora (na końcu linii z nazwą sekwencji/chromosomu). Wartość kursora powiększona o 1 odpowiada lokalizacji pierwszego nukleotydu tej sekwencji
								$line_len = 0; # Zmienna $line_len może zawierać dane ze wcześniejszych pętli, dlatego przed realizacją dalszej części kodu jest ona zerowana
								$alert = 0; # Zmienna $alert może zawierać dane ze wcześniejszych pętli, dlatego przed realizacją dalszej części kodu jest ona zerowana
							}
							else # Jeżeli linia nie rozpoczyna się od znaku zachęty '>', to ...
							{
								if ($line_len == 0) # ... sprawdza, czy zmienna $line_len jest równa 0 (pierwsza analizowana linia sekwencji). Jeżeli tak, to ...
								{
									$line_len = (split("", $_)); # ... do zmiennej $line_len przypisywana jest długość pierwszej linii z sekwencją nukleotydową
									$last_line_len = $line_len; # Zmiennej $last_line_len przypisywana jest wartość zmiennej $line_len (zmienne $last_line_len oraz $line_len mogą być sobie równe, jeżeli w referencji jest bardzo krótka sekwencja, lub linie z sekwencjami są bardzo długie)
									$chrom_len += $line_len; # Wartość zmiennej $chrom_len jest zwiększana o długość aktualnie analizowanej linii sekwencji
									print $Ofh "$line_len\t"; # Do pliku indeksu jest drukowana długość linii sekwencji nukleotydowej
								}
								else # Jeżeli zmienna $line_len jest różna od zera (jest to równoznaczne z tym, że są analizowane kolejne linie z sekwencjami nukleotydowymi), to ...
								{
									if ((split("", $_)) > 0) # ... sprawdza, czy dana linia nie jest pusta (czy są zapisane w niej jakiekolwiek znaki). Jeżeli tak, to ...
									{
										if ($line_len != (split("", $_))) # ... sprawdza, czy liczba znaków tej linii odpowiada liczbie znaków pierwszej analizowanej linii danej sekwencji/chromosomu. Jeżeli jest różna, to ...
										{
											if ($alert == 1) # ... sprawdza, czy zmienna $alert ma już przyporządkowaną wartość 1 (błąd odpowiadający różnej długości linii). Jeżeli tak, to ...
											{
												unlink $prefix . ".index"; # ... usuwany jest plik indeksu, następnie ...
												
												@genome_error = (2,"$chrom_ID,$.");
												last L;
											}
											$alert = 1; # Jeżeli zmienna $alert nie miała przyporządkowanej wcześniej wartości 1, to jest ona teraz jej nadawana
										}
										else
										{
											if ($alert == 1) # ... sprawdza, czy zmienna $alert ma już przyporządkowaną wartość 1 (błąd odpowiadający różnej długości linii). Jeżeli tak, to ...
											{
												unlink $prefix . ".index"; # ... usuwany jest plik indeksu, następnie ...
												
												@genome_error = (2,"$chrom_ID,$.");
												last L;
											}
										}
										$chrom_len += (split("", $_)); # Jeżeli wykonane wcześniej testy długości analizowanej linii wypadły korzystnie, to wartość zmiennej $chrom_len jest zwiększana o długość aktualnie analizowanej linii sekwencji oraz ...
										$last_line_len = (split("", $_)); # ... zmiennej $last_line_len przypisywana jest długość aktualnie analizowanej linii
									}
								}
							}
						}
						if ($chrom_len != 0) {print $Ofh "$last_line_len\t$chrom_len\n"; $chrom_len = 0}
					###########################
					close $fh;
					close $Ofh;
					
				} else {$genome_exists = 0}

				if ($genome_exists == 1 and $genome_error[0] == 0) { $reference_analysis_results[0] = 1 } # Jeżeli plik referencji istnieje ($genome_exists == 1) oraz nie ma żadnych błędów ($genome_error[0] == 0), to drukowany jest 'OK'
				elsif ($genome_exists == 0) { $reference_analysis_results[0] = 2 } # Jeżeli plik referencji nie istnieje ($genome_exists == 0), to drukowany jest 'fail'
				elsif ($genome_exists == 1 and $genome_error[0] == 1) { @reference_analysis_results = (3,$genome_error[1]) } # Jeżeli plik referencji istnieje ($genome_exists == 1), ale jest obecny w nim błąd ($genome_error[0] == 1), to drukowany jest 'fail' oraz nazwa zduplikowanej sekwencji ($genome_error[1])
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
		#	my %enzymes_for_analysis;

			if (defined $enzyme_file_name)
			{
				my $line = 0;
				my $enzymes_err = 0;
				my $enzymes_exists = 1;
				###################
				open $fh, '<', $enzyme_file_name;
					my @id_vs_companyName = (); # <ID_firmy>,<nazwa_firmy>
					
					
					while (<$fh>)
					{
						chomp $_;
						
						
						my @data = split(" ", $_);
						if (defined $data[0])
						{
							if ($date_check == 0) # Sprawdza datę bazy danych
							{
								for (my $i = 0; $i < scalar(@data); $i++)
								{
									my $check = 0;
									foreach my $month (qw/Jan Feb Mar Apr May Jun Jul Aug Sept Oct Nov Dec/)
									{
										if ($month eq $data[$i])
										{
											$check = 1;
											last;
										}
									}
									
									if ($check == 1)
									{
										$enzymes_db{date} = "$data[$i] $data[$i + 1] $data[$i + 2]";
										$date_check = 1;
									}
								}
							}
							elsif ($data[0] =~ /^[A-Z]$/ and $date_check == 1) # Zapisuje dane o producentach
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
								#$enzymes_db{companies}{$ID} = $company_name;
								$enzymes_db{companies} = join("\t",@id_vs_companyName);
							}
							elsif ($data[0] =~ /[A-Z][a-z]{2}/ and $date_check == 1)
							{
								$data[0] =~ s/;//;
								$data[2] =~ s/[_']//g;
								my @tmp = (); # <cut_site>,<sequence>,<overhang>,<isoschisomers>,<company_ID>
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
									#	$data[6] =~ s/>//;
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
								
								$enzymes_db{$data[0]} = join("\t",@tmp);
							#	print "$enzymes_db{$data[0]}\n";
							#	$enzymes_db{$data[0]}{cut_site} = $data[1];
							#	$enzymes_db{$data[0]}{sequence} = $data[2];
							#	$enzymes_db{$data[0]}{overhang} = $data[3];
								
								
								$data[2] = uc $data[2];
								push @enzymes_tmp, join("\t", $data[2],$data[0]);
							}				
						}
					}
				close $fh;
				

				
				if ($date_check == 0) { $enzymes_err = 1 }
				
				if ($enzymes_err == 0)
				{
					foreach my $enzyme (@enzymes_tmp)
					{
						push @enzymes_uniq_seq, ( split("\t", $enzyme) )[0]; # Tablica zawiera zestaw unikalnych sekwencji rozpoznawanych przez enzymy - do wykorzystania jako klucze %enzymes_for_analysis
					}
					@enzymes_uniq_seq = uniq(@enzymes_uniq_seq);
					
					foreach my $enzyme (@enzymes_tmp)
					{
						push @allEnzymesNames, split(",", ( split("\t", $enzyme) )[1] );  # Tablica zawiera zestaw unikalnych nazw enzymów - do wykorzystania jako klucze %enzymes_db
					}
					@allEnzymesNames = uniq(@allEnzymesNames);
					
					foreach my $enzyme_seq (@enzymes_uniq_seq) # Dla każdej unikalnej sekwencji rozpoznawanej przez enzym ...
					{
						my @data = ();
						my $seq_len = 0;
						for ( my $y = 0; $y < scalar(@enzymes_tmp); $y++ ) # Dla każdego rekordu enzymu pobranego z pliku gcg ...
						{
							if ( $enzyme_seq eq ( split("\t",$enzymes_tmp[$y]) )[0] ) # Jeżeli sekwencja unikalna jest identyczna z sekwencją z rekordu, to ...
							{
								push @data, ( split("\t",$enzymes_tmp[$y]) )[1]; # ... do tablicy @data dodawana jest nazwa enzymu z pasującego rekordu
								$seq_len = scalar( split("", ( split("\t",$enzymes_tmp[$y]) )[0] ) ); # Do zmiennej $seq_len zapisywana jest długość sekwencji z rekordu
							}
						}
						
						$enzymes_for_analysis{join(",",@data)} = join("\t",$enzyme_seq,$seq_len); # $enzymes_for_analysis{enz1,enz2,enz3,...} = ATGC\t4 (długość sekwencji)

					}
			
				}
				###################
				if ($enzymes_exists == 1 and $enzymes_err == 0) { @allEnzymesNames = sort{$a cmp $b} @allEnzymesNames; $enzyme_analysis_results[0] = 1 } # Jeżeli plik Enzymes istnieje ($enzymes_exists == 1) oraz nie ma żadnych błędów ($enzymes_err == 0), to drukowany jest 'OK'
				elsif ($enzymes_exists == 1 and $enzymes_err == 1) { $enzyme_analysis_results[0] = 3 } # Jeżeli plik Enzymes istnieje ($enzymes_exists == 1), ale jest obecny w nim błąd ($enzymes_err == 1), to drukowany jest 'fail' oraz jego lokalizacja (nr linii)
				
				
				
				
				
				print "Checking enzymes completed\n";

			}
			
			$jobID = 0;
		}
		elsif ($jobID == 3)
		{
			if (defined $vcf_file_name)
			{
				$line = 0;
				my $vcf_NoOfLines;
				my $vcf_error = 0; # Jeśli zmienna jest ustawiona na 1, to oznacza, że jest obecny błąd w pliku
				my $vcf_exists = 0;
				@linie = ();
				my $linie_ile = 0;
				
				if (-e $vcf_file_name) # Jeżeli plik VCF istnieje ...
				{
					$vcf_exists = 1;
					open my $check_vcf, '<', "$vcf_file_name";
						my $checkLine = 0;
						

						if ($vcf_exists == 1) # Jeżeli plik VCF istnieje, to ...
						{
							L: while (my $check_bier = <$check_vcf>)
							{
								chomp $check_bier;
								my @data = split("\t", $check_bier);
								$line++; # numeruje aktualnie sprawdzane linie w pliku
								if ($data[0] =~ /^#CHROM/ and $checkLine == 0) # Jeżeli plik VCF ma nazwy kolumn, to ...
								{
									for (my $i = 4; $i < scalar(@data); $i++) # ... do tablicy @linie zapisywane są nazwy analizowanych obiektów
									{
										push @linie, $data[$i];
									}
								}
								elsif ($checkLine == 0) {$vcf_error = 1; last L} # Jeżeli plik VCF nie ma nazw kolumn, to zmiennej $vcf_error przypisywana jest wartość 1, a następnie indukowane jest wyjście z pętli
								
								$linie_ile = @linie; # Przypisanie zmiennej $linie_ile liczby analizowanych obiektów

								if ($checkLine == 0) {$checkLine = 1}; # Zmienna $checkLine informuje czy aktualnie analizowany jest pierwszy wiersz zawierający nazwy kolumn ($checkLine == 0). Dla wszystkich pozostałych wierszy zmienna $checkLine ma przypisaną wartość 0

							}

							$vcf_NoOfLines = $line;
						}
					close $check_vcf;
				} else { $vcf_exists = 0 }
				$vcf_NoOfLines--;
				if ($vcf_exists == 1 and $vcf_error == 0) { %vcf_analysis_results = (err_code => 1, NoOfIndv => $linie_ile, NoOfSNPs => $vcf_NoOfLines) } # Jeżeli plik VCF istnieje ($vcf_exists == 1) oraz nie ma żadnych błędów ($vcf_error == 0), to drukowany jest 'OK'
				elsif ($vcf_exists == 0) { $vcf_analysis_results{err_code} = 2 } # Jeżeli plik VCF nie istnieje ($vcf_exists == 0), to drukowany jest 'fail'
				elsif ($vcf_exists == 1 and $vcf_error == 1) { $vcf_analysis_results{err_code} = 3 } # Jeżeli plik VCF istnieje ($vcf_exists == 1), ale nie ma nazw kolumn ($vcf_error == 1), to drukowany jest 'fail'
				elsif ($vcf_exists == 1 and $vcf_error == 2) { %vcf_analysis_results = (err_code => 4, errLine => $line) } # Jeżeli plik VCF istnieje ($vcf_exists == 1), ale ma co najmniej 1 niezdefiniowany genotyp "-" ($vcf_error == 2), to drukowany jest 'fail'

			}
			$jobID = 0;
		}
		elsif ($jobID == 4)
		{
			print "vcf2snps start\n";
			@samtools_error_code = (vcf2snps($snps_seq_len,$reference_file_name));
			
			$jobID = 0;
		}
		elsif ($jobID == 5)
		{
			print "snp2caps start\n";
			@snp2caps_results = ( snp2caps($output_seq_len,$reference_file_name) );
			
			$jobID = 0;
		}
		elsif ($jobID == 6)
		{
			print "singleCutSite start\n";
			@singleCutSite_results = ( singleCutSite() );
			
			$jobID = 0;
		}
		elsif ($jobID == 7)
		{
			print "Start simplifying VCF file\n";
			
			my $begin = 0;
			my @output;
			my @output_std = qw(0 1 3 4);
			my $vcf_error = 0;
			my $raw_vcf_file_name_tmp = $raw_vcf_file_name;
			$raw_vcf_file_name_tmp =~ s/.*[\\\/]//g;
			$raw_vcf_file_name_tmp =~ s/\..+$//g;
			$raw_vcf_file_name_tmp = $raw_vcf_file_name_tmp . ".svcf";
			
			open my $Ofh, '>>', $working_dir . $raw_vcf_file_name_tmp;
			open my $fh, '<', $raw_vcf_file_name;
				while (<$fh>)
				{
					chomp $_;
					my @data = split("\t", $_);
					if ($data[0] =~ /^#/ and $begin == 0)
					{
						if ($data[0] =~ /^#CHROM/ and $begin == 0)
						{
							$begin = 1;
						
							for (my $i = 9; $i < scalar(@data); $i++ )
							{
								push @output_std, $i;
							}
							
							print $Ofh join ("\t", @data[@output_std]);
							print $Ofh "\n";
						}						
					}
					elsif ($begin == 1)
					{
						for (my $i = 0; $i < 3; $i++ )
						{
							if (! print $Ofh "$data[ $output_std[$i] ]\t")
							{
								print "Error in writing to file $raw_vcf_file_name_tmp\n";
								exit;
							}
						}
						
						my @alt = split(",", $data[ $output_std[3] ]);
						print $Ofh "[";
						print $Ofh join(",", @alt);
						print $Ofh "]";
						
						for (my $i = 4; $i < scalar(@output_std); $i++ )
						{
							my @output = split(":", $data[ $output_std[$i] ]);	
							print $Ofh "\t$output[0]";
						}
						print $Ofh "\n";
					}
					elsif ($begin == 0)
					{
						$vcf_error = 1;
						last;
					}
				}
			close $fh;
			
			if ($vcf_error == 0) { %raw_vcf_analysis_results = (err_code => 1) }
			elsif ($vcf_error == 1) { %raw_vcf_analysis_results = (err_code => 2) }
			
			$jobID = 0;
		}
		sleep 1;
	}
	
	sub vcf2snps
	{
		print "vcf2snps start confirmed\n";
		print "\$working_dir: $working_dir\n";
		$samtools_error_code[1] = 0;
		my $genotypesLength = $_[0];
		my $error = 0;
	#	my $line = 0;
		$line_vcf = 0;
		@sequencesNotPresentInRef = ();
		$snpsNo = 0;
		@markersOnTheEdge = (); # Tablica zawierająca nazwy markerów, które są bliżej któregokolwiek końca niż wartość zawarta w polu 'DNA sequence length ...'
		
		my $prefix = $_[1];
		$prefix =~ s/\.[A-Za-z]+$//; # usunięcie suffixu '.fasta' lub '.fa' z nazwy pliku
		my %chrom_offset = (); # hash w którym będzie zapisywany ID chromosomu oraz pozycja kursora (offset)
		open $fh, '<', $prefix . ".index";
			while (<$fh>)
			{
				chomp $_;
				my @data = split("\t", $_);
				$data[0] =~ s/>//;
				$data[0] =~ s/ .+//g;
				$chrom_offset{$data[0]} = [($data[1],$data[2],$data[3],$data[4])];
			}
		close $fh;
		my $chrom_offset = \%chrom_offset;
		
		
		if (-e "$working_dir" . "snps.txt") {unlink "$working_dir" . "snps.txt"}

		open my $fh, '<', "$vcf_file_name";
		open my $Ofh, '>>', "$working_dir" . "snps.txt";

		while (my $bier = <$fh>)
		{
			$line_vcf++;

			chomp $bier;
			my @input = (split("\t", $bier)); # Zapis do tablicy @input wartości poszczególnych kolumn aktualnie analizowanego wiersza
			if ($input[0] eq "#CHROM") {next} # Jeżeli analizowana linia zawiera tekst '#CHROM', to jest ona pomijana
			my $RefSnpLen = scalar(split("", $input[2])); # Zapis do zmiennej długości sekwencji allelu SNPu z referencji
			my $ref_snp = $input[2]; # scalar zawierający SNP w sekwencji referencyjnej
			my $AltSnp = $input[3]; # Zapis do zmiennej $AltSnp wartości obecnej w 4 kolumnie zawierającej sekwencje wszystkiech alternatywnych alleli SNPu zamknięte w nawiasie kwadratowycm
			$AltSnp =~ s/[\[\]]//g; # Usunięcie ze zmiennej $AltSnp (zawierającej sekwencje wszystkiech alternatywnych alleli SNPu) nawiasów kwadratowych
			my @AltSnp = split(",", $AltSnp); # Zapis do tablicy @AltSnp sekwencji wszystkich alternatywnych alleli SNPu
			
			
			my $from = $input[1] - $genotypesLength; # Zapis do zmiennej $from numeru nukleotydu od którego będzie rozpoczynała się ekstrachowana sekwencja
			my $snp_min = $input[1] - 1; # Do zmiennej $snp_min zapisywana jest lokalizacja nukleotydu znajdującego się tuż przed SNP'em
			my $to = $input[1] + $genotypesLength + $RefSnpLen - 1; # Zapis do zmiennej $to numeru nukleotydu na którym będzie się kończyła ekstrachowana sekwencja. Zmienna $RefSnpLen zawiera informację o długości SNP'u referencji
			my $RefSnp_plus = $input[1] + $RefSnpLen; # Do zmiennej $RefSnp_plus zapisywana jest lokalizacja nukleotydu znajdującego się tuż za SNP'em. Zmienna $RefSnpLen zawiera informację o długości SNP'u referencji
			
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
			
			my $AltSnpNo = @AltSnp; # Zapis do zmiennej $AltSnpNo liczby alternatywnych alleli SNPu
			
			my @AltSnpLen = (); # tablica zawiera długości wszystkich alternatywnych wersji SNPów
			for (my $i = 0; $i < $AltSnpNo; $i++) # Każdy allelu SNPu ...
			{
				push @AltSnpLen, scalar(split("", $AltSnp[$i])); # ... jest dodawany do tablicy @AltSnpLen
			}

			
			

			
			
			
			
			my @genotypes = (); # Tablica zawierająca kody genotypów
			
			for (my $i = 4; $i < scalar(@input); $i++) # Dla każdej rośliny ...
			{
				if ( ( split("/", $input[$i]) )[0] eq "\." )
				{
					push @genotypes, $input[$i];
				}
				else
				{
					if ( ( split("/", $input[$i]) )[0] <= ( split("/", $input[$i]) )[1] ) # ... sprawdza, czy wartości liczbowe alleli w genotypie są ustawione w kolejności wzrastającej, tj. jeżeli np. genotyp jest równy 1/2, to ...
					{
						push @genotypes, $input[$i]; # ... do tablicy @geontypes dodawany jest dany genotyp
					}
					else # Jeżeli wartości liczbowe alleli w genotypie są ustawione w kolejności maljeącej, tj. jeżeli np. genotyp jest równy 2/1, to ...
					{
						push @genotypes, ( split("/", $input[$i]) )[1] . "/" . ( split("/", $input[$i]) )[0]; # ... do tablicy @geontypes dodawany jest dany genotyp o odwrotnej kolejności alleli, tj. 2/1 (przed) -> 1/2 (po)
					}
				}
				
				
			}
			

			my @samtoolsOut = samtools($reference_file_name,"$input[0]:$from-$snp_min",$chrom_offset); # pobiera całą sekwencję rejonu SNPu (z jego 'lewej' strony) z referencji


				
				
				
			print $Ofh ">$input[0]:$input[1]\n"; # Zapis identyfikatora SNPu do pliku wyjściowego 'out.txt'
			print $Ofh join("\t",@genotypes)."\n"; # Zapis genotypu każdej linii do pliku wyjściowego 'out.txt'
			
			my $extracted_sequence = $samtoolsOut[0];
		
		

			print $Ofh "$ref_snp,$RefSnpLen,$extracted_sequence$ref_snp"; # po wydrukowaniu 'lewej' części sekwencji SNPu drukuje sam SNP;
			
			@samtoolsOut = samtools($reference_file_name,"$input[0]:$RefSnp_plus-$to",$chrom_offset); # pobiera całą sekwencję rejonu SNPu (z jego 'prawej' strony) z referencji
			
			$extracted_sequence = $samtoolsOut[0];

			print $Ofh "$extracted_sequence";

			for (my $i = 0; $i < $AltSnpNo; $i++) # Dla każdego allelu SNPu ...
			{
				my $AltSnp_plus = $input[1] + $AltSnpLen[$i]; # ... do zmiennej $AltSnp_plus zapisywana jest lokalizacja nukleotydu znajdującego się tuż za SNP'em. Zmienna $AltSnpLen[$i] zawiera informację o długości danego allelu SNP'u
				$to = $input[1] + $genotypesLength + $AltSnpLen[$i] - 1; # Zapis do zmiennej $to numeru nukleotydu na którym będzie się kończyła ekstrachowana sekwencja. Zmienna $AltSnpLen[$i] zawiera informację o długości danego allelu SNP'u
				
				@samtoolsOut = samtools($reference_file_name,"$input[0]:$from-$snp_min",$chrom_offset); # pobiera całą sekwencję rejonu SNPu (z jego 'lewej' strony) z referencji
				
				$extracted_sequence = $samtoolsOut[0];

				print $Ofh "\t$AltSnp[$i],$AltSnpLen[$i],$extracted_sequence$AltSnp[$i]"; # po wydrukowaniu 'lewej' części sekwencji SNPu drukuje sam SNP;
				
				@samtoolsOut = samtools($reference_file_name,"$input[0]:$AltSnp_plus-$to",$chrom_offset); # pobiera całą sekwencję rejonu SNPu (z jego 'prawej' strony) z referencji
			
				$extracted_sequence = $samtoolsOut[0];

				print $Ofh "$extracted_sequence";
			}
			print $Ofh "\n";
			
			$snpsNo++;

		}
		print "vcf2snps finished\n";
		return ("",1);
	}
	
	sub snp2caps # Funkcja indentyfikująca SNP'y w sekwencjach rozpoznawanych przez enzymy restrykcyjne
	{
		my $seq_len = $_[0];
		my $seqName;
		my @genotypes;
		my %genotypes;
		my @seq;
		$stop = 0;
		@markersOnTheEdge = ();
		
		$numberOfSNPs = 0;
		my $linie_ile = scalar(@linie);
		print "Analysing ...\n\n";
		if (-e "$working_dir" . "out.txt") { unlink "$working_dir" . "out.txt" }
		
		$enzyme_zip = 0;

		my $prefix = $_[1];
		$prefix =~ s/\.[A-Za-z]+$//; # usunięcie suffixu '.fasta' lub '.fa' z nazwy pliku
		my %chrom_offset = (); # hash w którym będzie zapisywany ID chromosomu oraz pozycja kursora (offset)
		open $fh, '<', $prefix . ".index"; # otwarcie pliku indeksu
			while (<$fh>)
			{
				chomp $_;
				my @data = split("\t", $_);
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

		### UWAGA - rozważyć zastąpienie poniższej pętli zmienną $vcf_analysis_results{NoOfSNPs}
		##10 Zapis do zmiennej $fancy_numberOfSNPs liczby wszystkich SNP'ów obecnych w pliku 'snps.txt'
		open my $SNP, '<', "$working_dir" . "snps.txt"; # Otwarcie pliku 'snps.txt' (plik wsadowy powstały po przekształceniu pliku VCF)
			my $fancy_numberOfSNPs = 0; # Zmienna przechowująca liczbę wszystkich SNP'ów obecnych w pliku 'snps.txt'
			
			while (my $fancy_snp = <$SNP>)
			{
				chomp $fancy_snp;
				if ($fancy_snp =~ /^>/) # Jeżeli linia rozpoczyna się od znaku zachęty '>', to ...
				{
					$fancy_numberOfSNPs++; # ... wartość zmiennej $fancy_numberOfSNPs jest zwiększana o 1
				}
			}
			
		close $SNP;
		##10

		$fancy_SNPNo = 0; # Zmienna $fancy_SNPNo przechowuje nr aktualnie badanego SNP'u

		
		open $SNP, '<', "$working_dir" . "snps.txt"; # Otwarcie pliku 'snps.txt' (plik wsadowy powstały po przekształceniu pliku VCF)
		while (my $bierSNP = <$SNP>)
		{
			if ($stop == 1) { return ($numberOfSNPs,2) }
			
			chomp $bierSNP;

			my $IndvSeq = []; # Tablica zawierająca w każdej pozycji tablice z indywidualnymi nukleotydami fragmentu DNA z daną wersją alleliczną analizowanego SNP'u
			my @IndvSeq_snpLen_snp = (); # Tablica zawierająca długości każdej wersji allelicznej analizowanego SNP'u
			my @IndvSeq_snp = (); # Tablica zawierająca sekwencję każdej wersji allelicznej analizowanego SNP'u
			
			if ($bierSNP =~ />/) # Jeżeli linia rozpoczyna się od znaku zachęty '>', to ...
			{
				$seqName = $bierSNP; # ... zmiennej $seqName przyporządkowana jest wartość zmiennej $bierSNP przechowująca nazwę SNPu oraz ...
				$fancy_SNPNo++; # ... wartość zmiennej $fancy_SNPNo (nr aktualnie badanego SNP'u) jest zwiększana o 1
			}
			elsif ($bierSNP =~ /^[0-9\.]/) # Jeżeli linia rozpoczyna się od liczby, to ...
			{
				@genotypes = split("\t",$bierSNP); # ... do tablicy @genotypes zapisywane są kody liczbowe poszczególnych genotypów
				
				for (my $i = 0; $i < scalar(@linie); $i++) # Nazwom poszczególnych linii przypisywana jest wartość genotypu danego SNP'u
				{
					$genotypes{$linie[$i]} = $genotypes[$i];
				}
			}
			else # Jeżeli linia rozpoczyna się od sekencji SNP'u, to ...
			{
				@seq = split("\t",$bierSNP); # Zapis do tablicy @seq sekwencji SNP'u, jego długości oraz fragmentu DNA zawierającego dany SNP dla referencji oraz alternatywnych wersji SNP'u (fragment DNA zaczerpnięty z referencji, zmieniony jedynie nukleotyd w miejscu występowania SNP'u)
				my $SNP_alleles_No = @seq; # Zapis do zmiennej $SNP_alleles_No liczby alleli dango SNP'u
				
			#	my @enz = keys %enzymes_db; # Zapis do tablicy @enz nazw enzymów (dana pozycja może zawierać nazwę jednego lub większej liczby enzymów rozpoznających daną unikalną sekwencję DNA)
			#	@enz = grep { $_ ne "date" } @enz;
			#	@enz = grep { $_ ne "companies" } @enz;
				
				my $enz_No = @selected_enz_names_tmp; # Zapis do zmiennej $enz_No liczby unikalnych sekwencji rozpoznawanych przez enzymy
			#	print "\$enz_No: $enz_No\n";
			#	print join("\n", @selected_enz_names_tmp);
				
				
				
				for (my $i = 0; $i < $SNP_alleles_No; $i++) # Dla każdego allelu danego SNP'u ...
				{
					my @input = (split(",",$seq[$i])); # ... do tablicy @input zapisywane są następujące 3 dane: sekwencja danej wersji allelicznej SNP'u, długość SNP'u (danego allelu), fragment DNA ze SNP'em
					push @IndvSeq_snpLen_snp, $input[1]; # Do tablicy @IndvSeq_snpLen_snp zapisywane są długości danej wersji allelu danego SNP'u
					push @IndvSeq_snp, $input[0]; # Do tablicy @IndvSeq_snp zapisywana jest sekwencja danej wersji allelu danego SNP'u
					my @IndvSeq_tmp = split("",$input[2]); # Do tablicy @IndvSeq_tmp są zapisywane pojedyncze nukleotydy sekwencji danej wersji allelu danego SNP'u, a następnie ...
					push @$IndvSeq, [@IndvSeq_tmp]; # ... tablica ta jest zapisywana do tablicy zbiorczej @$IndvSeq
				}
				
				
				for (my $i = 0; $i < $enz_No; $i++) # Dla każdej unikalnej sekwencji rozpoznawanej przez enzymy restrykcyjne ...
				{			
					my $enzLength = split("", ( split("\t",$enzymes_db{$selected_enz_names_tmp[$i]}) )[1] ); # ... do zmiennej $enzLength zapisywana jest długość sekwencji rozpoznawanej przez enzym/enzymy
					#print "\$enzLength: $enzLength\n";
					$enzLength = $enzLength - 0; # ?

					my $cc = 0; # Zmienna $cc przechowująca kolejny nr wyekstrachowanego fragmentu DNA
					my @partSeq = ([],[]); # Tablica dwuwymiarowa przechowująca wyekstrachowane fragmenty DNA dla każdej wersji allelicznej badanego SNP'u. $partSeq[ID wersji allelicznej SNP'u][ID fragmentu DNA]
					my @partSeq_size = (); # ?
					my $numberOfMatches = 0; # ?
					my $numberOfMatches_all = 0; # ?
					my %matches1; # Hash przechowujący liczbę rozpoznawanych sekwencji (value) zawierających analizowany SNP dla każdej wersji allelicznej (key) analizowanego SNP'u
					my @partSeq1 = ""; # ?
					
					my ($to, $count) = ($enzLength,0); # Zmiennej $to przypisywana jest wartość zmiennej $enzLength zawierająca informacje o długości unikalnej sekwencji rozpoznawanej przez enzym. Zmiennej $count (?) przypisywana jest wartość 0.
					
					for (my $seqID = 0; $seqID < $SNP_alleles_No; $seqID++) # Dla każdej wersji allelicznej analizowanego SNP'u ...
					{
						$cc = 0; # ... resetowana (zerowana) jest zmienna $cc przechowująca kolejny nr wyekstrachowanego fragmentu DNA
						my $IndvSeq_Len = (scalar(@{$IndvSeq->[$seqID]})); # Przypisanie zmiennej $IndvSeq_Len długości fragmentu DNA zawierającego daną wersję alleliczną. Zmienna $seqID zawiera ID każdej wersji (0 -> referencja, 1 -> ...)
						my $from = int( ( ( $IndvSeq_Len - $IndvSeq_snpLen_snp[$seqID] ) + 1 ) / 2 ) - $enzLength; # Przypisanie zmiennej $from wartości, która w dalszej części programu jest wykorzystywana do obliczenia nr nukleotydu od którego mają być 'wycinane' fragmenty DNA o długości odpowiadającej długości sekwencji rozpoznawanej przez enzym i zawierających dany SNP
						

						
						##11 Ekstrakcja i zapis fragmentów DNA (do tablicy dwuwymiarowej $partSeq[ID danej wersji allelicznej][kolejny nr fragmentu]) o długości odpowiadającej długości sekwencji rozpoznawanej przez enzym i zawierających dany SNP
						for (my $c = 1; $c < $enzLength + ($IndvSeq_snpLen_snp[$seqID]); $c++) # Pętla ta określa ile fragmentów sekwencji ma zostać wyekstrachowanych (długość rozpoznawanej sekwencji + długość SNP'u)
						{
							my @part = (); # Tablica przechowująca pojedyncze nukleotydy wyekstrachowanego fragmentu DNA
						
							for (my $ii = $from + $c; $ii < $from + $enzLength + $c; $ii++)
							{

								my $input = ${@$IndvSeq[$seqID]}[$ii]; # Wyłuskiwanie pojedynczych nukleotydów (poczynając od wybranego początkowego nukleotydu) danego fragmentu DNA zawierającego daną wersję alleliczną badanego SNP'u ($seqID) do tymczasowej zmiennej $input
								push @part, $input; # Umieszczanie wyłuskanego nukleotydu w tablicy @part
								
							}
							
							$partSeq[$seqID][$cc] = join("",@part); # Łączenie pojedynczych nukleotydów w tablicy @part, a następnie zapisywanie ich do tablicy dwuwymiarowej $partSeq ($partSeq[ID wersji allelicznej SNP'u][ID fragmentu DNA])
							$cc++;
												
						}
						

						##11
					}

					##12 Zapis do tablicy @partSeq_size liczby wyekstrachowanych fragmentów DNA dla wszystkich wersji allelicznych danego SNP'u
					for (my $i = 0; $i < $SNP_alleles_No; $i++) # Dla każdej wersji allelicznej analizowanego SNP'u ...
					{
						my $input = @{$partSeq[$i]}; # ... do zmiennej $input przypisywana jest liczba wyekstrachowanych fragmentów DNA, a następnie ...
						@partSeq_size = (@partSeq_size, $input); # ... wartość tej zmiennej jest dodawana do tablicy @partSeq_size
					}
					##12
					my $enzyme_recogn_seq = ( split("\t",$enzymes_db{$selected_enz_names_tmp[$i]}) )[1];
					$enzyme_recogn_seq = uc $enzyme_recogn_seq;
					
					my $enzREGEX = enzREGEX( $enzyme_recogn_seq ); # przypisanie zmiennej $enzREGEX wyniku działania funkcji enzREGEX(<sekwencja DNA>), która zwraca sekwencję DNA w formie REGEX
					
					my $regex_inv = regex_inv($enzREGEX);
					if ($regex_inv eq $enzREGEX) # Sprawdzenie czy rozpoznawana przez testowany enzym sekwencja jest palindromem
					{
						##13 Określenie liczby wyekstrachowanych fragmentów DNA dla każdej wersji allelicznej analizowanego SNP'u, które są zrozpoznawane przez dany enzym restrykcyjny. Wartości te są zapisywane do hashu %matches1
						for (my $c = 0; $c < $SNP_alleles_No; $c++) # Dla każdej wersji allelicznej analizowanego SNP'u ...
						{ 
							for (my $i = 0; $i < $partSeq_size[$c]; $i++) # ... dla każdego wyekstrachowanego fragmentu DNA danej wersji allelicznej analizowanego SNP'u ...
							{	
								my $partSeq_upper = uc $partSeq[$c][$i];
								if ($partSeq_upper =~ /$enzREGEX/) # ... porównywana jest zgodność sekwencji danego fragmentu DNA z sekencją (w formie REGEX) rozpoznawaną przez enzym. Jeżeli sekwencje są identyczne, to ...
								{
									#	print "Sek: $c, len:  $partSeq[$c][$i] == $enzREGEX\n"; ####################################
									#	print "$enzyme_recogn_seq, $partSeq[$c][$i] == $enzREGEX\n";
									
									$numberOfMatches++; # ... wartość zmiennej $numberOfMatches (liczba sekwencji/miejsc ze SNP'em rozpoznawanych przez enzym dla danej wersji allelicznej analizowanego SNP'u) powiększana jest o 1
									$numberOfMatches_all++; # ... wartość zmiennej $numberOfMatches_all (liczba sekwencji/miejsc ze SNP'em rozpoznawanych przez enzym dla wszystkich wersji allelicznych analizowanego SNP'u) powiększana jest o 1
								}
		#							else
		#							{
		#								print "Sek: $c $partSeq[$c][$i] :: $enzREGEX\n"; ####################################
		#							}
							}
							$matches1{$c} = $numberOfMatches; # Do hashu %matches1 jest zapisywana liczba wyekstrachowanych fragmentów DNA ($numberOfMatches) dla każdej wersji allelicznej analizowanego SNP'u ($c), które są zrozpoznawane przez dany enzym restrykcyjny
							$numberOfMatches = 0; # Zmienna $numberOfMatches jest zerwowana przed kolejnym cyklem pętli
						}
						##13
					}
					else # Jeżeli rozpoznawana sekwencja nie jest palindromem, to testowany fragment DNA jest badany pod kątem występowania sekwencji komplementarnej do rozpoznawanej przez enzym
					{
						##13 Określenie liczby wyekstrachowanych fragmentów DNA dla każdej wersji allelicznej analizowanego SNP'u, które są zrozpoznawane przez dany enzym restrykcyjny. Wartości te są zapisywane do hashu %matches1
						for (my $c = 0; $c < $SNP_alleles_No; $c++) # Dla każdej wersji allelicznej analizowanego SNP'u ...
						{ 
							for (my $i = 0; $i < $partSeq_size[$c]; $i++) # ... dla każdego wyekstrachowanego fragmentu DNA danej wersji allelicznej analizowanego SNP'u ...
							{	
								my $partSeq_upper = uc $partSeq[$c][$i];
								if ($partSeq_upper =~ /$enzREGEX/ or $partSeq_upper =~ /$regex_inv/) # ... porównywana jest zgodność sekwencji danego fragmentu DNA z sekencją (w formie REGEX) rozpoznawaną przez enzym. Jeżeli sekwencje są identyczne, to ...
								{
									#	print "Sek: $c, len:  $partSeq[$c][$i] == $enzREGEX\n"; ####################################
									#	print "$enzyme_recogn_seq, $partSeq[$c][$i] == $enzREGEX\n";
									
									$numberOfMatches++; # ... wartość zmiennej $numberOfMatches (liczba sekwencji/miejsc ze SNP'em rozpoznawanych przez enzym dla danej wersji allelicznej analizowanego SNP'u) powiększana jest o 1
									$numberOfMatches_all++; # ... wartość zmiennej $numberOfMatches_all (liczba sekwencji/miejsc ze SNP'em rozpoznawanych przez enzym dla wszystkich wersji allelicznych analizowanego SNP'u) powiększana jest o 1
								}
		#							else
		#							{
		#								print "Sek: $c $partSeq[$c][$i] :: $enzREGEX\n"; ####################################
		#							}
							}
							$matches1{$c} = $numberOfMatches; # Do hashu %matches1 jest zapisywana liczba wyekstrachowanych fragmentów DNA ($numberOfMatches) dla każdej wersji allelicznej analizowanego SNP'u ($c), które są zrozpoznawane przez dany enzym restrykcyjny
							$numberOfMatches = 0; # Zmienna $numberOfMatches jest zerwowana przed kolejnym cyklem pętli
						}
						##13
					}
					
					


					
					##14 Określenie liczby genotypów/roślin, które są trawione przez enzym
					my $DigestedGenotypes_No_homo = 0; # Zmienna $DigestedGenotypes_No_homo przechowuje liczbę genotypów, dla których obydwa allele są trawione przez aktualnie badany enzym
					my $DigestedGenotypes_No_het = 0; # Zmienna $DigestedGenotypes_No_het przechowuje liczbę genotypów, dla których tylko jeden z dwóch alleli jest trawiony przez aktualnie badany enzym
					my $nullGenotypes = 0;
					for (my $i = 0; $i < $linie_ile; $i++) # Dla każdej analizowanej rośliny ...
					{
						my $genotypp = $genotypes{$linie[$i]}; # ... do zmiennej $genotypp jest zapisywany kod jej genotypu (0/0, 0/1, 1/2, itp.). Następnie ...
						my @genotypp_indv = split("/", $genotypp); # Do tablicy @genotypp_indv zapisywane są kody poszczególnych alleli SNPu (0 - ref, 1 - alternatywny, 2 - kolejny alternatywny, itd.)
						
						if ($genotypp_indv[0] ne "\.")
						{
							if ($genotypp_indv[0] == $genotypp_indv[1]) # Jeżeli obydwa allele SNPu obecne w danym genotypie są identyczne, to ...
							{
								if ($matches1{$genotypp_indv[0]} > 0) # ... sprawdza, ile wyekstrachowanych fragmentów (dla danego allelu SNPu) jest identyczna z sekwencją rozpoznawaną przez aktualnie badany enzym (kod genotypu odpowiada identyfikatorowi danej wersji allelicznej badanego SNP'u, np. genotyp 0 -> referencja, genotyp 1 -> pierwsza alternatywna wersja alleliczna, 2 -> ...)
								{
									$DigestedGenotypes_No_homo++; # Jeżeli warunek zostanie spełniony, to wartość zmiennej $DigestedGenotypes_No_homo jest zwiększana o 1 (dana roślina o danym genotypie ma allel w obrębie którego znajduje się miejsce rozpoznawane przez aktualnie badany enzym)
									#print "plant:$linie[$i]\tsingle:$genotypp\tmatches:$matches1{$genotypp}\n";
								}
							}
							elsif ($genotypp_indv[0] != $genotypp_indv[1]) # Jeżeli obydwa allele SNPu obecne w danym genotypie nie są identyczne, to ...
							{
								my $check_digest = 0;
								my $check_NO_digest = 0;
								foreach my $genot (@genotypp_indv) # ... dla każdego genotypu składowego ...
								{
									if ($matches1{$genot} > 0) # ... sprawdza, czy dla danego genotypu składowego (danej wersji allelicznej) stwierdzono obecność jakiejkolwiek wyekstrachowanego fragmentu DNA, który był identyczny z sekwencją rozpoznwaną przez enzym
									{
										$check_digest++;
									}
									else
									{
										$check_NO_digest++;
									}
								}
								
								if ($check_digest > 0 and $check_NO_digest == 0) # Jeżeli obydwa allele genotypu heterozygotycznego ulegają trawieniu, to ...
								{
									$DigestedGenotypes_No_homo++; #  ... wartość zmiennej $DigestedGenotypes_No_homo jest zwiększana o 1 (obydwa allele genotypu heterozygotycznego danej rośliny mają miejsca rozpoznawane przez aktualnie badany enzym)
								}
								elsif ($check_digest > 0 and $check_NO_digest > 0) # Jeżeli tylko jeden allel genotypu heterozygotycznego ulega trawieniyu, to ...
								{
									$DigestedGenotypes_No_het++; #  ... wartość zmiennej $DigestedGenotypes_No_het jest zwiększana o 1 (tylko jeden allel genotypu heterozygotycznego danej rośliny ma miejsce rozpoznawane przez aktualnie badany enzym -> na żelu będą 2 fragmenty)
								}
		#						print "plant:$linie[$i]\tduo:$genotypp\n";
							}
						}
						else
						{
							$nullGenotypes++;
						}
					}
	#				print $DigestedGenotypes_No_homo . "\n";
	#				print $DigestedGenotypes_No_het . "\n";
	#				print join(" ",keys %genotypes);
	#				exit;
					##14

					##15 Zapis do pliku 'out.txt' fragmentu sekwencji z daną wersją alleliczną analizowanego SNP'u (dla każdej rośliny)
					my @genotypes_uniq = uniq(@genotypes); # Do tablicy @genotypes_uniq zapisywane są wszystkie rodzaje genotypów występujące u badanych obiektów
					my @genotypes_uniq_tmp = grep { $_ ne "\.\/\." } @genotypes_uniq; # Zapisuje do tablicy @genotypes_uniq_tmp listę unikalnych genotypów za wyjątkiem './.'. Do dalszych analiz będą brane jedynie te markery, które mają co najmniej 2 różne genotypy (nie uwzględniając './.')
					

					if ($polymorphicSNPsOnly == 1 and $numberOfMatches_all > 0 and $DigestedGenotypes_No_homo < ( $linie_ile - $nullGenotypes ) and scalar(@genotypes_uniq_tmp) > 1 ) # (Aktywna opcja 'Mine CAPS from polymorphic ...') Jeżeli pośród analizowanych roślin są genotypy trawione i nietrawione, to wynik jest zapisywany do pliku
					{ 
						$numberOfSNPs++;

						my ($chrom,$snp) = split(":",$seqName); # Ekstrakcja z nazwy analizowanego markera/SNP'u nazwy sekwencji/chromosomu na którym jest zlokalizowany SNP ($chrom) oraz nr nukleotydu w obrębie tej sekwencji/chromosomu w którym jest zlokalizowany dany SNP lub od którego rozpoczyna się sekwencja indelu ($snp)
						my $snp_min = $snp - 1; # Do zmiennej $snp_min zapisywana jest lokalizacja nukleotydu znajdującego się tuż przed SNP'em
						$chrom =~ s/>//; # Usunięcie z nazwy sekwencji/chromosomu znaku zachęty '>'
						my $taker_from = $snp - $seq_len; # Zapis do zmiennej $taker_from numeru nukleotydu od którego będzie rozpoczynała się ekstrachowana sekwencja
						my $snp_plus = $snp + $IndvSeq_snpLen_snp[0]; # Do zmiennej $snp_plus zapisywana jest lokalizacja nukleotydu znajdującego się tuż za SNP'em. Zmienna $IndvSeq_snpLen_snp[0] zawiera informację o długości SNP'u referencji
						my $taker_to = $snp + $seq_len + $IndvSeq_snpLen_snp[0] - 1; # Zapis do zmiennej $taker_to numeru nukleotydu na którym będzie się kończyła ekstrachowana sekwencja. Zmienna $IndvSeq_snpLen_snp[0] zawiera informację o długości SNP'u referencji
						
						if ( $taker_from < 0 or $taker_to > $chrom_offset->{$chrom}[3] )
						{
							push @markersOnTheEdge, $seqName;
							$numberOfSNPs--;
							next;
						}
						
						open my $fh, '>>', "$working_dir" . "out.txt"; # Otwarcie do zapisu pliku wynikowego 'out.txt'
							print $fh "$seqName\n";
							my $enz_seq = ( split("\t",$enzymes_db{$selected_enz_names_tmp[$i]}) )[1];
							my $tmp = $selected_enz_names_mix{$selected_enz_names_tmp[$i]};
							print $fh "$tmp\t$enz_seq\n";
							
							
							#my $taker_to = $snp + $seq_len - ($IndvSeq_snpLen_snp[0] - 1); # Zapis do zmiennej $taker_to numeru nukleotydu na którym będzie się kończyła ekstrachowana sekwencja. Zmienna $IndvSeq_snpLen_snp[0] zawiera informację o długości SNP'u referencji
							
	#						print "\n";
	#						print "\$taker_to = $snp + $seq_len - ($IndvSeq_snpLen_snp[0] - 1) = $taker_to";
	#						exit;
							
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
								
								my @taker = samtools($reference_file_name,"$chrom:$taker_from-$snp_min",$chrom_offset); # Ekstrakcja 'lewej' części sekwencji (jej ostatni nukleotyd znajduje się tuż przez SNP'em referencji)
								$taker[0] = lc $taker[0]; # Zamiana dużych liter wyekstrachowanej sekwencji na małe
								print $fh $taker[0] . $IndvSeq_snp[$i]; # Zapis do pliku 'out.txt' 'lewej' części sekwencji ($taker[0]) oraz sekwencji samego SNP'u ($IndvSeq_snp[0])
								
								@taker = samtools($reference_file_name,"$chrom:$snp_plus-$taker_to",$chrom_offset); # Ekstrakcja 'prawej' części sekwencji (jej pierwszy nukleotyd znajduje się tuż za SNP'em referencji)
								$taker[0] = lc $taker[0]; # Zamiana dużych liter wyekstrachowanej sekwencji na małe
								print $fh "$taker[0]\n"; # Zapis do pliku 'out.txt' 'prawej' części sekwencji ($taker[0])
							}
							
							#my @genotypes_uniq = uniq(@genotypes); # Do tablicy @genotypes_uniq zapisywane są wszystkie rodzaje genotypów występujące u badanych obiektów

							foreach my $genotype_uniq ( sort{$b cmp $a} @genotypes_uniq ) # Dla każdego rodzaju genotypu ...
							{
								print $fh "$genotype_uniq\t"; # ... zapis do pliku wynikowego 'out.txt' nazwy genotypu, następnie ...
								print $fh "["; # ... zapis do pliku wynikowego 'out.txt' nawiasu '['
								
								my $c = 0; # Zmienna określająca z który allel genotypu jest aktualnie analizowany (0 - allel pierwszy, 1 - allel drugi)
								foreach my $indv_allele ( split("/", $genotype_uniq) ) # Dla każdego allelu wchodzącego w skład danego genotypu ...
								{
									if ($indv_allele eq "\.")
									{
										print $fh "?";
									}
									else
									{
										if ($matches1{$indv_allele} == 0) # ... jeżeli nie znajduje się w miejscu rozpoznawanym przez enzym, to ...
										{
											print $fh "-"; # ... zapisywany do pliku wynikowego 'out.txt' jest znak '-'
										}
										else
										{
											print $fh "+"; # A jeżeli dany allel znajduje się w miejscu rozpoznawanym przez enzym, to durkowany jest znak '+'
										}
									}
									
									if ($c == 0) # Jeżeli aktualnie analizowany jest pierwszy allel genotypu, to ...
									{
										print $fh "/"; # ... drukowany jest znak '/'
									}
									else
									{
										print $fh "]"; # ... a jeśli nie, to drukowany jest znak ']'
										
										foreach my $indv_name (sort{$a cmp $b} keys %genotypes) # Dla posortowanej nazwy każdego obiektu ...
										{										
											if ($genotypes{$indv_name} eq $genotype_uniq) # ... jeżeli genotyp odpowiadający obiektowi o danej nazwie jest identyczny z aktualnie tetsowanym genotypem, to ...
											{
												print $fh "\t$indv_name"; # ... do pliku wynikowego drukowana jest nazwa tego obiektu
											}
										}
									}
									
									$c++; # Zwiększenie wartości zmiennej $c o 1 (sygnał, że w kolejnej rundzie będzie analizowany drugi allel)
								}
								print $fh "\n"; # Po przeanalizowaniu danego rodzaju genotypu drukowany jest znak nowej linii '\n'
							}
							print $fh "\n";
						close $fh;
					}
					elsif ($polymorphicSNPsOnly == 0 and $numberOfMatches_all > 0 and ( $DigestedGenotypes_No_homo > 0 or $DigestedGenotypes_No_het > 0 ) ) # (Nieaktywna opcja 'Mine CAPS from polymorphic ...') Jeżeli genotyp którejkolwiek rośliny jest trawiony, to wynik jest zapisywany do pliku
					{
						$numberOfSNPs++;

						my ($chrom,$snp) = split(":",$seqName); # Ekstrakcja z nazwy analizowanego markera/SNP'u nazwy sekwencji/chromosomu na którym jest zlokalizowany SNP ($chrom) oraz nr nukleotydu w obrębie tej sekwencji/chromosomu w którym jest zlokalizowany dany SNP lub od którego rozpoczyna się sekwencja indelu ($snp)
						my $snp_min = $snp - 1; # Do zmiennej $snp_min zapisywana jest lokalizacja nukleotydu znajdującego się tuż przed SNP'em
						$chrom =~ s/>//; # Usunięcie z nazwy sekwencji/chromosomu znaku zachęty '>'
						my $taker_from = $snp - $seq_len; # Zapis do zmiennej $taker_from numeru nukleotydu od którego będzie rozpoczynała się ekstrachowana sekwencja
						my $snp_plus = $snp + $IndvSeq_snpLen_snp[0]; # Do zmiennej $snp_plus zapisywana jest lokalizacja nukleotydu znajdującego się tuż za SNP'em. Zmienna $IndvSeq_snpLen_snp[0] zawiera informację o długości SNP'u referencji
						my $taker_to = $snp + $seq_len + $IndvSeq_snpLen_snp[0] - 1; # Zapis do zmiennej $taker_to numeru nukleotydu na którym będzie się kończyła ekstrachowana sekwencja. Zmienna $IndvSeq_snpLen_snp[0] zawiera informację o długości SNP'u referencji
						
						if ( $taker_from < 0 or $taker_to > $chrom_offset->{$chrom}[3] )
						{
							push @markersOnTheEdge, $seqName;
							$numberOfSNPs--;
							next;
						}
						
						open my $fh, '>>', "$working_dir" . "out.txt"; # Otwarcie do zapisu pliku wynikowego 'out.txt'
							print $fh "$seqName\n";
							my $enz_seq = ( split("\t",$enzymes_db{$selected_enz_names_tmp[$i]}) )[1];
							my $tmp = $selected_enz_names_mix{$selected_enz_names_tmp[$i]};
							print $fh "$tmp\t$enz_seq\n";
							
							
							#my $taker_to = $snp + $seq_len - ($IndvSeq_snpLen_snp[0] - 1); # Zapis do zmiennej $taker_to numeru nukleotydu na którym będzie się kończyła ekstrachowana sekwencja. Zmienna $IndvSeq_snpLen_snp[0] zawiera informację o długości SNP'u referencji
							
	#						print "\n";
	#						print "\$taker_to = $snp + $seq_len - ($IndvSeq_snpLen_snp[0] - 1) = $taker_to";
	#						exit;
							
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
								
								my @taker = samtools($reference_file_name,"$chrom:$taker_from-$snp_min",$chrom_offset); # Ekstrakcja 'lewej' części sekwencji (jej ostatni nukleotyd znajduje się tuż przez SNP'em referencji)
								$taker[0] = lc $taker[0]; # Zamiana dużych liter wyekstrachowanej sekwencji na małe
								print $fh $taker[0] . $IndvSeq_snp[$i]; # Zapis do pliku 'out.txt' 'lewej' części sekwencji ($taker[0]) oraz sekwencji samego SNP'u ($IndvSeq_snp[0])
								
								@taker = samtools($reference_file_name,"$chrom:$snp_plus-$taker_to",$chrom_offset); # Ekstrakcja 'prawej' części sekwencji (jej pierwszy nukleotyd znajduje się tuż za SNP'em referencji)
								$taker[0] = lc $taker[0]; # Zamiana dużych liter wyekstrachowanej sekwencji na małe
								print $fh "$taker[0]\n"; # Zapis do pliku 'out.txt' 'prawej' części sekwencji ($taker[0])
							}
							
							#my @genotypes_uniq = uniq(@genotypes); # Do tablicy @genotypes_uniq zapisywane są wszystkie rodzaje genotypów występujące u badanych obiektów

							foreach my $genotype_uniq ( sort{$b cmp $a} @genotypes_uniq ) # Dla każdego rodzaju genotypu ...
							{
								print $fh "$genotype_uniq\t"; # ... zapis do pliku wynikowego 'out.txt' nazwy genotypu, następnie ...
								print $fh "["; # ... zapis do pliku wynikowego 'out.txt' nawiasu '['
								
								my $c = 0; # Zmienna określająca z który allel genotypu jest aktualnie analizowany (0 - allel pierwszy, 1 - allel drugi)
								foreach my $indv_allele ( split("/", $genotype_uniq) ) # Dla każdego allelu wchodzącego w skład danego genotypu ...
								{
									if ($indv_allele eq "\.")
									{
										print $fh "?";
									}
									else
									{
										if ($matches1{$indv_allele} == 0) # ... jeżeli nie znajduje się w miejscu rozpoznawanym przez enzym, to ...
										{
											print $fh "-"; # ... zapisywany do pliku wynikowego 'out.txt' jest znak '-'
										}
										else
										{
											print $fh "+"; # A jeżeli dany allel znajduje się w miejscu rozpoznawanym przez enzym, to durkowany jest znak '+'
										}
									}
									
									if ($c == 0) # Jeżeli aktualnie analizowany jest pierwszy allel genotypu, to ...
									{
										print $fh "/"; # ... drukowany jest znak '/'
									}
									else
									{
										print $fh "]"; # ... a jeśli nie, to drukowany jest znak ']'
										
										foreach my $indv_name (sort{$a cmp $b} keys %genotypes) # Dla posortowanej nazwy każdego obiektu ...
										{										
											if ($genotypes{$indv_name} eq $genotype_uniq) # ... jeżeli genotyp odpowiadający obiektowi o danej nazwie jest identyczny z aktualnie tetsowanym genotypem, to ...
											{
												print $fh "\t$indv_name"; # ... do pliku wynikowego drukowana jest nazwa tego obiektu
											}
										}
									}
									
									$c++; # Zwiększenie wartości zmiennej $c o 1 (sygnał, że w kolejnej rundzie będzie analizowany drugi allel)
								}
								print $fh "\n"; # Po przeanalizowaniu danego rodzaju genotypu drukowany jest znak nowej linii '\n'
							}
							print $fh "\n";
						close $fh;
					}
					
					

					

					#parser($IndvSeq,$SNP_alleles_No,$enzymes{$selected_enz_names[$i]},$enzLength,$seqName,$selected_enz_names[$i],[%genotypes],[@IndvSeq_snpLen_snp],[@IndvSeq_snp]);
					# $IndvSeq - tabela z sekwencjami każdej wersji genotypu
					# $SNP_alleles_No - liczba genotypów
					# $enzymes{$selected_enz_names[$i]} - sekwencja rozpoznawana przez enzym
					# $enzLength - długość sekwencji rozpoznawanej przez enzym
					# $seqName - identyfikator SNPu
					# $selected_enz_names[$i] - nazwa enzymu
					# @genotypes - lista ID genotypów (0 - ref., [1-9] - altern.)
					# @IndvSeq_snpLen_snp - długość danego SNPu
					# @IndvSeq_snp - sekwencja SNPu
				}
			}
		}
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
		$numberOfSNPsBefore = 0;
		$numberOfSNPsAfter = 0;
		$stop = 0;
		
		if (-e "$working_dir" . "out_single.txt") { unlink "$working_dir" . "out_single.txt" }
		
		open my $Ofh, '>>', "$working_dir" . "out_single.txt";
		open my $fh, '<', "$working_dir" . "out.txt";
			while (my $bier = <$fh>)
			{
				if ($stop == 1) { return ("",2) }
				
				chomp $bier;
				if ($bier =~ /^>/) # Jeżeli linia w pliku wynikowym 'out.txt' rozpoczyna się od znaku zachęty '>', to ...
				{
					if ( scalar(keys %AllSeq) != 0 ) # Jeżeli hash %AllSeq zawiera jakiekolwiek dane ...
					{
						my $check = 0;
						foreach my $snp_allele_ID (keys %AllSeq) # Dla każdego allelu SNPu ...
						{
							if ($AllSeq_NoOfRecognSeq{$snp_allele_ID} > 1) # ... sprawdza, czy liczba sekwencji rozpoznawanych przez enzym jest równa 1. Jeżeli tak, to ...
							{
								$check++;
							}
						}
						
						if ($check == 0)
						{
							foreach my $line (@output) # ... każda linia rekordu danego SNPu zapisaną w tablicy @output jest zapisywana do pliku wynikowego 'out_single.txt'
							{
								print $Ofh "$line\n";
							}
							
							print $Ofh "\n";
							$numberOfSNPsAfter++; # Wartość zmiennej $numberOfSNPsAfter (liczba SNPów z allelami zawierający) jest zwiększana o 1
						}
					}
					
					
					%AllSeq = ();
					%AllSeq_NoOfRecognSeq = ();
					@output = ();
					
					push @output, $bier; # ... do tablicy @output zapisywana jest cała linia oraz ...
					$numberOfSNPsBefore++; # ... wartość zmiennej $numberOfSNPsBefore zwiększana jest o 1
				}
				elsif ($bier =~ /^[A-Z]/ and scalar(split("\t", $bier)) == 2 ) # Jeżeli linia w pliku wynikowym 'out.txt' rozpoczyna się od dużej litery (linia z nazwami enzymów), to ...
				{
					@input = split("\t",$bier); # Do dablicy @input zapisywane są poszczególne pola danej linii
					$EnzSeq = $input[1];
					push @output, $bier; # ... do tablicy @output zapisywana jest cała linia oraz ...
				}
				elsif ($bier =~ /^[ra]/ and scalar(split("\t", $bier)) == 5) # Jeżeli linia w pliku wynikowym 'out.txt' rozpoczyna się od 'ref', to ...
				{
					push @output, $bier; # ... do tablicy @output zapisywana jest cała linia oraz ...
					@input = split("\t",$bier); # Do dablicy @input zapisywane są poszczególne pola danej linii
					$preAllSeq = $input[4]; # Do zmiennej $preAllSeq zapisywana jest sekwencja nukleotydowa
					$preAllSeq = uc $preAllSeq; # Zamiana dużych liter na małe w zmiennej $preAllSeq
					my $snp_ID = $input[1]; # Do zmiennej $snp_ID zapisywany jest identyfikator danego allelu SNPu
					$snp_ID =~ s/[\[\]]//g; # Usunięcie z wartości zmiennej $snp_ID nawiasów kwadratowych
					$AllSeq{$snp_ID} = $preAllSeq; # Zapis do hashu %AllSeq sekwencji ($preAllSeq) zawierającej każdy allel SNPu
					#push @AllSeq, $preAllSeq;
				}
				elsif ($bier =~ /^[0-9\.]/ and scalar(split("\t", $bier)) >= 3)
				{
					push @output, $bier; # ... do tablicy @output zapisywana jest cała linia
				}
				
				$EnzSeq = uc $EnzSeq;
				my $enzREGEX = enzREGEX($EnzSeq);
				

				@EnzSeqSingleNucl = split("",$EnzSeq);
				$EnzSeqSingleNuclSize = @EnzSeqSingleNucl; 
				
				foreach my $snp_allele_ID (keys %AllSeq) # Dla każdego allelu SNPu ...
				{
					my $NoOfRecognSeq = 0;
					my @Seq_singleNucl = split("", $AllSeq{$snp_allele_ID});
					my $Seq_singleNucl_size = @Seq_singleNucl;
					
					## Pętla dzieli analizowaną sekwencję z danym allelem SNPu na którsze fragmenty o długości równej długości sekwencji rozpoznawanej przez enzym, a następnie porównuje te fragmenty z sekwencją ropoznawaną przez enzym.
					my $regex_inv = regex_inv($enzREGEX);
					if ($regex_inv eq $enzREGEX) # Sprawdzenie czy rozpoznawana przez testowany enzym sekwencja jest palindromem
					{
						for (my $i = 0; $i < $Seq_singleNucl_size - $EnzSeqSingleNuclSize; $i++)
						{
							my @seqPart = @Seq_singleNucl[$i..( $i + ($EnzSeqSingleNuclSize - 1) )]; # Ekstrakcja fragmentu o długości równej długości sekwencji rozpoznawanej przez enzym
							my $seqPart = join("",@seqPart);
							my $seqPart_upper = uc $seqPart;
							if (uc $seqPart_upper =~ /$enzREGEX/)
							{
								$NoOfRecognSeq++;
							}
						}
						##
					}
					else
					{
						for (my $i = 0; $i < $Seq_singleNucl_size - $EnzSeqSingleNuclSize; $i++)
						{
							my @seqPart = @Seq_singleNucl[$i..( $i + ($EnzSeqSingleNuclSize - 1) )]; # Ekstrakcja fragmentu o długości równej długości sekwencji rozpoznawanej przez enzym
							my $seqPart = join("",@seqPart);
							my $seqPart_upper = uc $seqPart;
							if (uc $seqPart_upper =~ /$enzREGEX/ or uc $seqPart_upper =~ /$regex_inv/)
							{
								$NoOfRecognSeq++;
							}
						}
						##
					}
					
					
					
					
					$AllSeq_NoOfRecognSeq{$snp_allele_ID} = $NoOfRecognSeq; # Zapis do hashu %AllSeq_NoOfRecognSeq liczby rozpoznawanych sekwencji dla danego allelu SNPu
				}
			}
		close $fh;
		
		
		
		## Fragment kodu odpowiedzialny za sprawdzenie, czy dany CAPS (z pliku 'out.txt') zawiera sekwencję z jakimkolwiek allelem SNP z co najwyżej 1 miejscem trawienia. Jeżeli prawda, to dany CAPS jest drukowany do pliku 'out_single.txt'
		my $check = 0;
		foreach my $snp_allele_ID (keys %AllSeq) # Dla każdego allelu SNPu ...
		{
			if ($AllSeq_NoOfRecognSeq{$snp_allele_ID} > 1) # ... sprawdza, czy liczba sekwencji rozpoznawanych przez enzym jest równa 1. Jeżeli tak, to ...
			{
				$check++;
			}
		}
		
		if ($check == 0)
		{
			foreach my $line (@output) # ... każda linia rekordu danego SNPu zapisaną w tablicy @output jest zapisywana do pliku wynikowego 'out_single.txt'
			{
				print $Ofh "$line\n";
			}
			
			print $Ofh "\n";
			$numberOfSNPsAfter++; # Wartość zmiennej $numberOfSNPsAfter (liczba SNPów z allelami zawierający) jest zwiększana o 1
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
	
	sub samtools
	{
		my ($fh,$Ofh);
		my $error_code = 1;
		
		
		my $input = $_[0]; # nazwa pliku fasta
		my $chrom_and_coord = $_[1]; # <ID chromosomu>:<nr nukleotydu początkowego>-<nr nukleotydu końcowego>
		my $chrom = (split(":", $chrom_and_coord))[0]; # zapis ID chromosomu
		my $coord = (split(":", $chrom_and_coord))[1]; # zapis koordynatów ekstrachowanej sekwencji
		my $from_nucl = (split("-", $coord))[0]; # nr nukleotydu początkowego
		my $to_nucl = (split("-", $coord))[1]; # nr nukleotydu końcowego
		my $chrom_offset_ = $_[2];
		
		if (!$chrom_offset_->{$chrom})
		{
			return ($chrom,3);
		}
		elsif ($from_nucl eq "" or $to_nucl > $chrom_offset_->{$chrom}[3])
		{
		#	print "$coord $to_nucl\n";
			return ($chrom,4);
		}
		

		open $fh, '<', $input or die "Cannot open the file $input\n";
			my $chrom_offset = $chrom_offset_->{$chrom}[0];
			my $line_len = $chrom_offset_->{$chrom}[1];
			my $last_line_len = $chrom_offset_->{$chrom}[2];
			my $chrom_len = $chrom_offset_->{$chrom}[3];

			if (!defined($chrom_offset))
			{
#				print "Error - there is no '$chrom' sequence in the '$input' file.\n";
#				print "Exiting ...\n";
				$error_code = 3;
				return ($chrom,$error_code);
			}
			
			my $from_offset = 0;
			my $offset = 0;
			
			if ($from_nucl <= ($chrom_len - $last_line_len + 1) or $chrom_len == $last_line_len)
			{
				$from_offset = $chrom_offset + $from_nucl - 1 + int($from_nucl / $line_len); # pozycja kursora od którego rozpocznie się czytanie sekwencji
				$offset = $to_nucl - $from_nucl + 1 + (int($to_nucl / $line_len) - int($from_nucl / $line_len));
				if (($from_nucl % $line_len) == 0) {$from_offset -= 1; $offset += 1} # parametr ustalony eksperymentalnie :/
			}
			else
			{
				$from_offset = $chrom_offset + $from_nucl + int($from_nucl / $line_len); # pozycja kursora od którego rozpocznie się czytanie sekwencji
				$offset = $to_nucl - $from_nucl + (int($to_nucl / $line_len) - int($from_nucl / $line_len));
				if (($from_nucl % $line_len) != 0) {$from_offset -= 1; $offset += 1} # parametr ustalony eksperymentalnie :/
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
}