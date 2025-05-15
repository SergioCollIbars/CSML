#!/usr/bin/perl
 
=pod

=begin HEADER

=head1 NAME

    hpf_eef_transform.pl - transforms (XML based) EEF into plain text formats

=head1 SYNOPSIS

    Usage: hpf_eef_transform.pl [--help] [--verbose] [--nocleanup] [--config <config file>] [--ouput <output dir>] <EEF file(s)>

=head1 DESCRIPTION

    See manual (GO-TN-HPF-GS-0192)

=head1 COPYRIGHT

    COPYRIGHT (c) 2008 - 2011 SRON (Netherlands Institute for Space Research)

    This is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License, version 2, as
    published by the Free Software Foundation.

    The software is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 51 Franklin Street, Fifth Floor, 
    Boston, MA  02110-1301, USA.

=head1 AUTHOR

    S. de Witte (hpf-2000@sron.nl)
    S. Fiorot (hpf-2000@sron.nl)
    R. van Hees (hpf-2000@sron.nl)

=head1 VERSION

    1.2.1 25-Feb-2011 Fixed several Window's incompatibilities
                      Fixed problems on small memory PC's with EEF (level 1b) products
    1.2.0 23-Dec-2010 Parse and transform Level 1b products (EGG_NOM_1b, SST_NOM_1b, STR_VC2/3_1b)
    1.1.6 19-Nov-2010 modify the match test so that no xml declaration is copied in the .tmp file, also 
                      if a word wrap (return) sign (from windows ?) is in the HDR file.
    1.1.5 12-Nov-2010 update Max_Recursion_Depth (hpf_eef_transform.cfg). Increase to 1999.
                      add a default path to the .cfg directory. If the --config not used, then the default .cfg
		      file is located in the same directory as the .pl file 
                      add a default path for the stylesheets if the stylesheets directory (<Location>) is not
		      specified in the .cfg file then the stylesheets are in the same directory as the .pl file.
    1.1.4 27-Sep-2010 update spatial GG format
    1.1.3 01-Jun-2010 update time GG format
    1.1.1 31-May-2010 bug fixed: output_dir parsing at line 149
    1.1   10-Jul-2008 introduced new XSLT modes when transforming the
                      last splitted filename (to avoid writing "end of
		      file" markers in the middle of an output file)
    1.0.1 07-Jul-2008 bugfix when running without arguments
    1.0   20-Jun-2008 first official release
    0.4   16-Jun-2008 decoding of "Base64" encoded binary data added
    0.3   13-Jun-2008 changed "splitted" tag into "Split"
    0.2   12-May-2008 extra tag within encoded products
    0.1    9-May-2008 initial release

=end HEADER
=cut

use strict;
use warnings;

use Getopt::Long;
use File::Basename;
use XML::Parser;
use XML::LibXML;
use XML::LibXSLT;
use MIME::Base64 qw( decode_base64 );

# set constant config file XML tags
# (note: one can use the unary plus operator '+'
# when using one of the constants as a hash key)
use constant {
    CONFIG_ROOT         => 'HPF_EEF_Transform_Configuration',
    CONFIG              => 'Configuration',
    PRODUCTS            => 'Products',
    PRODUCT             => 'Product',
    NAME                => 'Name',
    HAS_SUBPRODUCTS     => 'Has_Subproducts',
    EXTENSIONS          => 'Extensions',
    INPUT               => 'Input',
    OUTPUT              => 'Output',
    HEADER              => 'Header',
    DATA                => 'Data',
    TEMP                => 'Temp',
    MAX_RECURSION_DEPTH => 'Max_Recursion_Depth',
    INTEGER             => 'Integer',
    STYLESHEETS         => 'Stylesheets',
    LOCATION            => 'Location'
};
    
# declare global variables
our ($g_version, $g_program, $g_verbose, $g_cleanup);
our ($g_dir_separator, $g_single_dir_separator, $g_multiple_dir_separators);
our ($g_xmlfile, $g_datasets, $g_output_products );
our (%g_file_handle, %g_joint_datafile_info);

# set current version
unless (defined $g_version) {
    my $major = '1';
    my $minor = '2';
    my $built = '1';
    $g_version = "$major.$minor.$built";
}

# set name of current program
( $g_program = $0 ) =~ s/.*\///o;

# define the directory separator and regular expressions
# for a single and multiple directory separator(s)
$g_dir_separator           = '/';
$g_single_dir_separator    = qr/$g_dir_separator/;
$g_multiple_dir_separators = qr/$g_dir_separator{2,}/;


MAIN:
{

    # force a buffer flush after every print
    $| = 1;

    # define default configuration file
    my $default_config_extension = '.cfg';
    my ($program_filename_without_extension, undef, undef) = fileparse($g_program, qr/\.[^.]*/);
    my $program_path=dirname($0);
    my $default_config_file = $program_path .'/' . $program_filename_without_extension . $default_config_extension;

    # declare and empty local variable
    # for storing command line options
    my %opt = ();

    # get command line options and store in hash
    GetOptions(\%opt, "help|?", "verbose", "nocleanup", "config=s", "output=s");

    # print usage info to STDERR when asked for help
    usage($default_config_file) if (defined $opt{help});

    # set verbosity level
    $g_verbose = $opt{verbose} || 0;

    # define whether to delete temporary splitted files
    $g_cleanup = $opt{nocleanup} ? 0 : 1;

    # define the configuration file
    my $config_file = $opt{config} || $default_config_file;

    # define and check output directory (if specified),
    # make sure that it ends with a single dir separator
    my $output_dir = $opt{output} if (defined $opt{output});
    if (defined $output_dir) {
        $output_dir =~ s/$g_multiple_dir_separators/$g_dir_separator/g;
        $output_dir =~ s/${g_single_dir_separator}$//;
        die if (dir_is_not_ok($output_dir, 'output'));
        $output_dir =~ s/$/${g_dir_separator}/;
    }

    # the rest of the command line should be input filename(s)
    # (in EEF/XML format), at least one filenames should be given
    @ARGV >= 1 or usage($default_config_file);

    # read the configuration file
    my $config = parse_config_file($config_file);

    # check the filenames from the command line
    # and store them in a reference to a hash
    my ($checked_filenames, $header_found, $datafile_found);
    if ((defined $config) && (ref $config eq "HASH") && (exists $config->{+EXTENSIONS})) {
        ($checked_filenames, $header_found, $datafile_found) = check_command_line_filenames(\@ARGV, $config->{+EXTENSIONS}, $config->{+PRODUCTS});
    }
    
    #Check the path of the stylesheets
    unless ( defined $config->{+STYLESHEETS}->{+LOCATION} ) {
         $config->{+STYLESHEETS}->{+LOCATION} = $program_path;
    } 

    # define target size for splitted XML files
    my $percentage_of_memory = 7;
    my $target_size = determine_target_size($percentage_of_memory);

    # transform each file group from the command line
    my $process_result;
    foreach my $input_product_instance (keys %{$checked_filenames}) {
      
        if ( (($header_found->{$input_product_instance}) && ($datafile_found->{$input_product_instance})) || 
             (($datafile_found->{$input_product_instance})  && ($checked_filenames->{$input_product_instance}->{+DATA}->{ext} eq '.EEF')) ) {

	    # if an output directory is not specified, check
	    # whether we can write in the input file directories
	    unless (defined $output_dir) {

                # define file types using constants
                my @file_types;
                if (($checked_filenames->{$input_product_instance}->{+DATA}->{ext} eq '.EEF')) {
                   # here we just have 1 file with both header info and data_block info
                   @file_types = (DATA);
                } else {
                   @file_types = (DATA, HEADER);
                }
		# define directories to be checked
		# and store them as keys in a hash
		# overwriting identical directories
		my %dirs_to_be_checked;
		foreach my $file_type (@file_types) {
		    if (exists $checked_filenames->{$input_product_instance}->{$file_type}->{dir}) {
		        $dirs_to_be_checked{$checked_filenames->{$input_product_instance}->{$file_type}->{dir}} = 1;
		    }
		}

		# check directories and exit when they are not writable
	        foreach my $input_directory_used_as_output_directory (keys %dirs_to_be_checked) {
		    die if (dir_is_not_ok($input_directory_used_as_output_directory, 'output'));
		}

	    }

	    # everything ok, so process files for current product instance
	    my $max_recursion_depth = $config->{+MAX_RECURSION_DEPTH} || 1999;
	    $process_result = process_files({prod_inst           => $input_product_instance,
			                     filenames           => $checked_filenames->{$input_product_instance},
			                     extensions          => $config->{+EXTENSIONS},
			                     products            => $config->{+PRODUCTS},
			                     max_recursion_depth => $max_recursion_depth,
			                     output_dir          => $output_dir,
			                     target_size         => $target_size,
			                     stylesheets_dir     => $config->{+STYLESHEETS}->{+LOCATION}});

	} else {

	    # print missing file(s) error
	    my %error;
	    $error{type} = "ERROR";
            if ($header_found->{$input_product_instance}) {
	        $error{message} = DATA;
	    } elsif ($datafile_found->{$input_product_instance}) {
	        $error{message} = HEADER;
	    }
	    $error{message} .= "file(s) missing on command line for $input_product_instance";
	    if ((exists $config->{+EXTENSIONS}) && (exists $config->{+EXTENSIONS}->{+INPUT})) {
	        if (($header_found->{$input_product_instance}) && (defined $config->{+EXTENSIONS}->{+INPUT}->{+HEADER})) {
	            $error{message} .= $config->{+EXTENSIONS}->{+INPUT}->{+HEADER};
	        } elsif (($datafile_found->{$input_product_instance}) && (defined $config->{+EXTENSIONS}->{+INPUT}->{+DATA})) {
	            $error{message} .= $config->{+EXTENSIONS}->{+INPUT}->{+DATA};
	        }
	    }
            print_error(\%error);

	}

    }

    # exit gracefully
    if (defined $process_result) {
        print STDOUT "Result: $process_result\n" if ($g_verbose);
        exit 0;
    } else {
        print STDOUT "Result: errors encountered\n" if ($g_verbose);
        exit 1;
    }

}


sub parse_config_file {

    my $config_file = shift;

    # declare return variable
    my $configuration;

    # check whether there is a configuration file and whether it's OK
    return if (( ! defined $config_file) || (file_is_not_ok($config_file, 'configuration')));

    # define a new XML parser object
    my $config_parser_options = config_parser_options();
    my $config_parser         = XML::Parser->new(%{$config_parser_options}) if (defined $config_parser_options);

    # parse the configuration file
    $configuration = $config_parser->parsefile($config_file);

    # return configuration
    return $configuration;

}


sub config_parser_options {

    # define handlers
    my $handlers = {Start => \&config_parser_start_tag,
                    Char  => \&config_parser_contents,
                    Final => \&config_parser_results};
    
    return {Handlers => $handlers};

}


sub config_parser_start_tag {

    my ($config_parser, $element, %attributes) = @_;

    if ((defined $element) && ($config_parser->within_element(CONFIG_ROOT)) && ($config_parser->within_element(CONFIG))) {
        if (($config_parser->within_element(PRODUCTS)) && ($element eq PRODUCT)) {
	    if ((exists $attributes{+NAME}) && (exists $attributes{+HAS_SUBPRODUCTS})) {
                $config_parser->{+CONFIG}->{+PRODUCTS}->{$attributes{+NAME}} = $attributes{+HAS_SUBPRODUCTS};
	    } elsif (exists $attributes{+NAME}) {
                $config_parser->{+CONFIG}->{+PRODUCTS}->{$attributes{+NAME}} = 0;
            }
	}
    }

}


sub config_parser_contents {

    my ($config_parser, $contents) = @_;

    if ((defined $contents) && ($config_parser->within_element(CONFIG_ROOT)) && ($config_parser->within_element(CONFIG))) {
	if ($config_parser->within_element(MAX_RECURSION_DEPTH)) {
	    my $max_recursion_depth = trim($contents);
	    $max_recursion_depth = 1999 unless ($max_recursion_depth =~ /^\d+$/o);
	    $config_parser->{+CONFIG}->{+MAX_RECURSION_DEPTH} = $max_recursion_depth if ($config_parser->current_element eq INTEGER);
	} elsif ($config_parser->within_element(EXTENSIONS)) {
	    my $extension = trim($contents);
	    $extension = '.' . $extension unless ($extension =~ /^\./o);
	    if ($config_parser->within_element(INPUT)) {
	        $config_parser->{+CONFIG}->{+EXTENSIONS}->{+INPUT}->{+HEADER} = $extension if ($config_parser->current_element eq HEADER);
	        $config_parser->{+CONFIG}->{+EXTENSIONS}->{+INPUT}->{+DATA}   = $extension if ($config_parser->current_element eq DATA);
	    } elsif ($config_parser->within_element(OUTPUT)) {
	        $config_parser->{+CONFIG}->{+EXTENSIONS}->{+OUTPUT}->{+HEADER} = $extension if ($config_parser->current_element eq HEADER);
	        $config_parser->{+CONFIG}->{+EXTENSIONS}->{+OUTPUT}->{+DATA}   = $extension if ($config_parser->current_element eq DATA);
	        $config_parser->{+CONFIG}->{+EXTENSIONS}->{+OUTPUT}->{+TEMP}   = $extension if ($config_parser->current_element eq TEMP);
	    }
	} elsif ($config_parser->within_element(STYLESHEETS)) {
	    my $location = trim($contents);
	    $location .= $g_dir_separator unless ($location =~ /${g_single_dir_separator}$/);
	    $config_parser->{+CONFIG}->{+STYLESHEETS}->{+LOCATION} = $location if ($config_parser->current_element eq LOCATION);
	}
    }

}


sub config_parser_results {

    my ($config_parser) = @_;

    # declare return variable
    my %configuration;

    # store the configuration in a hash
    if (exists $config_parser->{+CONFIG}) {
        %configuration = %{$config_parser->{+CONFIG}};
    }

    # return a reference to the (sorted) splitted filenames
    return \%configuration;

}


sub determine_target_size {

    my $phys_mem_percentage = shift;

    # declare return variable
    my $target_size;
    
    # initialize default memory variables
    my $phys_mem_total = 128;
    my $phys_mem_unit  = 'Mb';

    # get memory info
    my %memory_info = get_mem_info();
    if ((exists $memory_info{MemTotal}) && ($memory_info{MemTotal} =~ /^[+-]?\d+\.?\d*$/o)) {
        if ((exists $memory_info{MemTotalUnit}) && ($memory_info{MemTotalUnit} =~ /^\s*[kmg]+b?\s*$/io)) {
            $phys_mem_total = $memory_info{MemTotal};
            $phys_mem_unit  = $memory_info{MemTotalUnit};
	}
    }

    # define part of physical memory
    my $phys_mem_partial = int(($phys_mem_total * $phys_mem_percentage) / 100);

    # define target size
    my %byte_factor = ( ' ' => 1, K => 1_024, M => 1_048_576, G => 1_073_741_824);
    if ($phys_mem_unit =~ /^\s*([kmg]+)b?\s*$/io) {
	my $unit = "\U$1";
        $target_size = $phys_mem_partial * $byte_factor{$unit};
    }

    # return target size
    return $target_size;

}


sub get_mem_info {

    # declare return vartiable
    my %mem;

    if (open(MEMINFO, "< /proc/meminfo")) {

        foreach (<MEMINFO>) {
            if (/^Mem:\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)/) {
                $mem{mem_total}   = $1;
                $mem{mem_used}    = $2;
                $mem{mem_free}    = $3;
                $mem{mem_shared}  = $4;
                $mem{mem_buffers} = $5;
                $mem{mem_cached}  = $6;
                next;
            } elsif (/^Swap:\s+(\S+)\s+(\S+)\s+(\S+)/) {
                $mem{swap_total} = $1;
                $mem{swap_used}  = $2;
                $mem{swap_free}  = $3;
                next;
            } elsif(/^(\S+):\s+(\S+) (\S+)/) {
                my $unit    = $1 . "Unit";
                $mem{$1}    = $2;
                $mem{$unit} = $3;
                next;
            }
        }

        close(MEMINFO);

        return %mem;

    } else {
        
	return;

    }

}


sub process_files {

    my ($arguments) = @_;

    my $product_instance    = $arguments->{prod_inst};
    my $input_files         = $arguments->{filenames};
    my $extension           = $arguments->{extensions};
    my $products            = $arguments->{products};
    my $max_recursion_depth = $arguments->{max_recursion_depth};
    my $ouput_dir           = $arguments->{output_dir};
    my $target_size         = $arguments->{target_size};
    my $stylesheets_dir     = $arguments->{stylesheets_dir};

    # declare return variable
    my $transform_results;

    # declare a hash for storing the full path
    # to the temporary output file names
    my %tmp_output_filename_full_path;
        
    # declare full path and input- and ouput-filenames variables
    my ($input_header_dir, $output_header_dir, $input_data_file_dir, $output_data_file_dir);
    my ($input_header_full_path, $input_data_file_full_path);

    # if Level 2 files (== .HDR and .DBL in the command line)
    if (defined $input_files->{+HEADER}->{ext}) {
	# define input and output directories
	$input_header_dir     = $input_files->{+HEADER}->{dir};
	$output_header_dir    = (defined $ouput_dir) ? $ouput_dir : $input_header_dir;
	$input_data_file_dir  = $input_files->{+DATA}->{dir};
	$output_data_file_dir = (defined $ouput_dir) ? $ouput_dir : $input_data_file_dir;

	# define input header filename
	$input_header_full_path = $input_header_dir . $g_dir_separator . $product_instance . $input_files->{+HEADER}->{ext};
	$input_header_full_path =~ s/$g_multiple_dir_separators/$g_dir_separator/g;

	# define input data filename
	$input_data_file_full_path = $input_data_file_dir . $g_dir_separator . $product_instance . $input_files->{+DATA}->{ext};
	$input_data_file_full_path =~ s/$g_multiple_dir_separators/$g_dir_separator/g;
	   
	# make sure we know where to read from and where to write to
	return unless ((defined $input_header_full_path) && (defined $output_header_dir));
	return unless ((defined $input_data_file_full_path) && (defined $output_data_file_dir));

	# if desired, print some extra information before processing
	if ($g_verbose) {
	    print STDOUT "INPUT HEADER    : $input_header_full_path\n";
	    print STDOUT "INPUT DATA FILE : $input_data_file_full_path\n";
	} 
	# check whether everyting is OK with the input files
	my $input_file_error = 0;
	$input_file_error = 1 if (file_is_not_ok($input_header_full_path, 'header'));
	$input_file_error = 1 if (file_is_not_ok($input_data_file_full_path, 'data'));
	return if ($input_file_error);
    } else {   
        # Level 1 b data : only 1 EEF file in the command line
        # define input and output directories
        $input_data_file_dir  = $input_files->{+DATA}->{dir};
        $output_data_file_dir = (defined $ouput_dir) ? $ouput_dir : $input_data_file_dir;

        # define input data filename
        $input_data_file_full_path = $input_data_file_dir . $g_dir_separator . $product_instance . $input_files->{+DATA}->{ext};
        $input_data_file_full_path =~ s/$g_multiple_dir_separators/$g_dir_separator/g;
        return unless ((defined $input_data_file_full_path) && (defined $output_data_file_dir));

        if ($g_verbose) {
	    print STDOUT "INPUT DATA FILE : $input_data_file_full_path\n";
        }

        # check whether everyting is OK with the input files
        my $input_file_error = 0;
        $input_file_error = 1 if (file_is_not_ok($input_data_file_full_path, 'data'));
        return if ($input_file_error);
    }
    
    # check the EEF file against the filename convention and determine its properties;
    # skip this EEF file when the product is not defined, in that case return undefined
    my ($mission, $file_class, $eef_product, undef, $validity_period, $file_version, undef) = parse_eef_filename($input_data_file_full_path);
    return unless ((defined $mission) && (defined $file_class) && (defined $eef_product) && (defined $validity_period) && (defined $file_version));

    # define splitting parameters for splitting the datafile
    my $split_level = get_split_level($eef_product, $products);
    my $number_of_digits_for_splitted_filenames = 4;
    my $splitted_filename_extension = $extension->{+OUTPUT}->{+TEMP} || '.tmp';
    my $splitted_filename_options = {digits => $number_of_digits_for_splitted_filenames,
                                     ext    => $splitted_filename_extension};
    
    # define the expected (after a successfull split) XML index file (with the splitted filenames within their (sub)product tags)
    my $expected_index_file_full_path = $output_data_file_dir . $g_dir_separator . $product_instance . '-' . '0' x $number_of_digits_for_splitted_filenames . $splitted_filename_extension;
    $expected_index_file_full_path =~ s/$g_multiple_dir_separators/$g_dir_separator/g;


    # split the EEF into HDR and DBL files
    if ($input_files->{+DATA}->{ext} eq '.EEF' ) {
	#define usual directories to split the hdr and dbl files      
	$output_data_file_dir = (defined $ouput_dir) ? $ouput_dir : $input_data_file_dir;
	$output_header_dir    = $output_data_file_dir;
	$input_header_dir     = $output_data_file_dir;
	$input_data_file_dir  = $output_data_file_dir;

	# output DBL file should be of nearly equal size as original	
	my $target_size_dbl = -s $input_data_file_full_path;

	my $split_result_eef = split_xml({data_file        => $input_data_file_full_path,
					  header           => $input_header_full_path,
					  output_dir       => $output_data_file_dir,
					  target_size      => $target_size_dbl,
					  level            => $split_level,
					  filename_options => $splitted_filename_options});
     

	# move the 2 splitted files to the so-called .HDR and .DBL.
	# They are located in the output directory
	my $splitted_hdr_file = $output_data_file_dir . $g_dir_separator . $product_instance . '-0001.tmp';
	my $splitted_dbl_file = $output_data_file_dir . $g_dir_separator . $product_instance . '-0002.tmp';
	my $new_hdr_file = $output_data_file_dir . $g_dir_separator . $product_instance . '.HDR';
	my $new_dbl_file = $output_data_file_dir . $g_dir_separator . $product_instance . '.DBL'; 
        my $new_file = $output_data_file_dir . $g_dir_separator . $product_instance;

	rename( $splitted_hdr_file, $new_hdr_file );
	rename( $splitted_dbl_file, $new_dbl_file );
	unlink( "$new_file*.tmp" );

	#define usual path to split the hdr and dbl files
	$input_header_full_path   = $new_hdr_file;
	$input_data_file_full_path= $new_dbl_file;
    }

    # split the datafile into smaller parts (so they will fit into memory for the XSL transformation)
    my $split_result = split_xml({data_file        => $input_data_file_full_path,
                                  header           => $input_header_full_path,
				  output_dir       => $output_data_file_dir,
				  target_size      => $target_size,
				  level            => $split_level,
				  filename_options => $splitted_filename_options});

    # check whether the expected XML index file was created
    # and, if so, assign it to the actual index file variable
    my $index_file_full_path;
    if ($split_result) {
        $index_file_full_path = $expected_index_file_full_path unless ((file_is_not_ok($expected_index_file_full_path, 'index')) || ( -z $expected_index_file_full_path));
    } else {
    	my %error;
	$error{type} = "ERROR";
	$error{message} = "splitting data file ($input_data_file_full_path) failed";
        print_error(\%error);
	return;
    }

    # parse the XML index file
    my $splitted_filenames;
    if (defined $index_file_full_path) {
        $splitted_filenames = parse_index_file($index_file_full_path);
    }
    
    # transform files using XSL stylesheets
    if ((defined $splitted_filenames) && (ref $splitted_filenames eq "HASH") && (exists $extension->{+OUTPUT})) {
        
	# cleanup of temporary index file
	if (($g_cleanup) && (-e $index_file_full_path)) {
	    unless (unlink $index_file_full_path) {
	        my $system_error = $!;
		my ($filename, $dir, $ext) = fileparse($index_file_full_path, qr/\.[^.]*/);
		my %error;
		$error{type}    = "WARNING";
		$error{message} = "cannot delete temporary index file";
		$error{message} .= " (" . $filename . $ext . ") ";
		$error{message} .= "from \"$dir\"";
		$error{message} .= ": $system_error" if (defined $system_error);
		print_error(\%error);
	    }
	}

	# define unique ID's for header and data block
	# (note: these are used within the XSL stylesheets as well)
	my $last_input_file_indicator = 'last_file';
	my $header                    = 'header';
	my $header_and_last_file      = $header . '_and_' . $last_input_file_indicator;
	my $data_block                = 'data_block';
	my $data_block_and_last_file  = $data_block . '_and_' . $last_input_file_indicator;

	# store mode id's in hash
	my %mode_ids;
	$mode_ids{'header'}               = $header;
	$mode_ids{'header_last_file'}     = $header_and_last_file;
	$mode_ids{'data_block'}           = $data_block;
	$mode_ids{'data_block_last_file'} = $data_block_and_last_file;
	
	# define stylesheets (case sensistive)
	my %stylesheet;
	my $stylesheet_common_name = 'eef2txt';
	my $stylesheet_extension   = '.xsl';

	# define header stylesheet
	if ($eef_product =~ /_(1)(b|B)$/ ) {
	    $stylesheet{$header} = $stylesheets_dir . $g_dir_separator . $stylesheet_common_name . '_' . $header . '_L1b' . $stylesheet_extension;
            $stylesheet{$header} =~ s/$g_multiple_dir_separators/$g_dir_separator/g;
        } else {    
	    $stylesheet{$header} = $stylesheets_dir . $g_dir_separator . $stylesheet_common_name . '_' . $header . $stylesheet_extension;
            $stylesheet{$header} =~ s/$g_multiple_dir_separators/$g_dir_separator/g;
	}
	
	
	# define data block stylesheet based on the EEF product
        $stylesheet{$data_block} = $stylesheets_dir . $g_dir_separator . $stylesheet_common_name . '_' . $data_block . '_' . $eef_product . $stylesheet_extension;
        $stylesheet{$data_block} =~ s/$g_multiple_dir_separators/$g_dir_separator/g;

        # transform files
        $transform_results = transform_files({filenames            => $splitted_filenames,
	                                      input_product_inst   => $product_instance,
	                                      stylesheet           => \%stylesheet,
					      mode_ids             => \%mode_ids,
					      max_recusrion_depth  => $max_recursion_depth,
					      output_ext           => $extension->{+OUTPUT},
					      output_header_dir    => $output_header_dir,
					      output_data_file_dir => $output_data_file_dir,
					      mission              => $mission,
					      file_class           => $file_class,
					      validity_period      => $validity_period,
					      file_version         => $file_version});
    }
    
 
    $index_file_full_path =~ /(.*)-(\d{4})\.tmp/;
    my $file_withoutext = $1;
        
    # delete (temporarly created) HDR and DBL in case of Level 1b 
    if ( ($g_cleanup) && ($index_file_full_path =~ /GO_(.*)_1b(.*)/) ) {
        if ($g_verbose) {
	    print STDOUT "Delete temporary files : $file_withoutext.HDR $file_withoutext.DBL\n";
        }
      	unless( unlink "$file_withoutext.HDR" ) {
	    my $system_error = $!;
	    my ($filename, $dir, $ext) = fileparse( "$file_withoutext.HDR", qr/\.[^.]*/);
	    my %error;
	    $error{type}    = "WARNING";
	    $error{message} = "cannot delete temporary EEF header-files";
	    $error{message} .= " (" . $filename . $ext . ") ";
	    $error{message} .= "from \"$dir\"";
	    $error{message} .= ": $system_error" if (defined $system_error);
	    print_error(\%error);
	}
      	unless( unlink "$file_withoutext.DBL" ) {
	    my $system_error = $!;
	    my ($filename, $dir, $ext) = fileparse( "$file_withoutext.DBL", qr/\.[^.]*/);
	    my %error;
	    $error{type}    = "WARNING";
	    $error{message} = "cannot delete temporary EEF data-file";
	    $error{message} .= " (" . $filename . $ext . ") ";
	    $error{message} .= "from \"$dir\"";
	    $error{message} .= ": $system_error" if (defined $system_error);
	    print_error(\%error);
	}
    }
    
    # Only some flags (ccd, ggt, iaq) are extracted from the Level 1b data : but
    # there is a tmp file  created for each flag : so we need to remove them at the end of the XML Parser
    if ( ($g_cleanup) ) {
	my @not_extracted_file = glob($file_withoutext.'*.tmp');
	foreach my $not_extracted_file_to_be_deleted (@not_extracted_file) {
	    if ($g_verbose) {
		print STDOUT "Delete intermediate files : $not_extracted_file_to_be_deleted\n";
	    }
	    unlink( "$not_extracted_file_to_be_deleted" );
	}
    }

    return $transform_results;

}


sub get_split_level {

    my $eef_product         = shift;
    my $configured_products = shift;

    # declare return variable
    my $split_level;

    # check whether the EEF product has subproducts
    my $eef_product_has_subproducts = 0;
    foreach my $configured_product (keys %{$configured_products}) {
        if ($eef_product =~ /$configured_product/i) {
	    $eef_product_has_subproducts = 1 if ($configured_products->{$configured_product});
	    last;
	}
    }

    # define split level
    if ($eef_product =~ /_(1)(b|B)$/ ) {
       $split_level = 2
    } else {
       $split_level = ($eef_product_has_subproducts) ? 4 : 3;
    }

    return $split_level;

}


sub split_xml {

    my ($arguments) = @_;

    my $data_file_full_path = $arguments->{data_file};
    my $header_full_path    = $arguments->{header};
    my $output_dir          = $arguments->{output_dir};
    my $target_size         = $arguments->{target_size};
    my $split_level         = $arguments->{level};
    my $filename_options    = $arguments->{filename_options};

    # declare return variable
    my $split_result;

    # make sure everything is defined
    return unless ( ((defined $data_file_full_path) && (defined $header_full_path)) ||
                    ((defined $data_file_full_path) && ($data_file_full_path =~ /.EEF$/ )) );
    return unless ((defined $output_dir) && (defined $target_size) && (defined $split_level) && (defined $filename_options));

    # split header (if not .EEF the first time we enter here) and data file names
    my ($header_name, $header_dir, $header_ext)          = fileparse($header_full_path, qr/\.[^.]*/) if ($data_file_full_path !~ /.EEF$/);
    my ($data_file_name, $data_file_dir, $data_file_ext) = fileparse($data_file_full_path, qr/\.[^.]*/);

    # define a new XML split object using the filename options
    my $splitted_xml = xml_split::splitted_xml::split_parser->new($filename_options);
    $splitted_xml->{current_size}   = 0;
    $splitted_xml->{target_size}    = $target_size;
    $splitted_xml->{handlers}       = {Start => \&split_parser_start_tag, End => \&split_parser_end_tag, Default => \&split_parser_default};
    $splitted_xml->{target_depth}   = $split_level;
    $splitted_xml->{seq_number}     = 0;
    $splitted_xml->{header_name}    = $header_name;
    $splitted_xml->{header_dir}     = $header_dir;
    $splitted_xml->{header_ext}     = $header_ext;
    $splitted_xml->{new_root}       = 'Earth_Explorer_File';
    $splitted_xml->{data_file_name} = $data_file_name;
    $splitted_xml->{data_file_dir}  = $data_file_dir;
    $splitted_xml->{data_file_ext}  = $data_file_ext;
    $splitted_xml->{output_dir}     = $output_dir;


    # define a new XML parser object with options
    my $split_parser_options = split_parser_options($splitted_xml);
    my $split_parser         = XML::Parser->new(%{$split_parser_options}) if (defined $split_parser_options);

    if (defined $split_parser) {
    
        # add a reference to the element stack array
        # (needed to keep track of open elements) and
        # the parser object to the XML split object 
        my @element_stack;
        $splitted_xml->{stack}  = \@element_stack;
        $splitted_xml->{parser} = $split_parser;

        # if desired, print some extra information during splitting
	my $start_info_message_length;
	if ($g_verbose) {
	    my $info_message = "Splitting input data file into smaller files (target size: $target_size bytes). This can take a while, please wait...";
	    $start_info_message_length = length($info_message);
            print STDOUT "$info_message\r";
	}
	
        # parse and split the data file
        $split_result = $split_parser->parsefile($data_file_full_path);
	
        # if desired, print some extra information when done splitting
	# (extra spaces at the end are needed to clear previous message)
	my $end_info_message_length;
	if ($g_verbose) {
            my $padding_character = ' ';
	    my $info_message = "Finished splitting input data file.";
	    $end_info_message_length = length($info_message);
	    if ((defined $start_info_message_length) && (defined $end_info_message_length) && ($start_info_message_length > $end_info_message_length)) {
	        $info_message .= ($padding_character x ($start_info_message_length - $end_info_message_length));
            }
	    print STDOUT "$info_message\n";
	}
    }
         
    # close the index filehandle (if still open)
    if (defined $splitted_xml->{index_fh}) {
        unless (close($splitted_xml->{index_fh})) {
            my $system_error = $!;
	    my %error;
    	    $error{type}    = "WARNING";
	    $error{message} = "cannot close index file";
	    $error{message} .= ": $system_error" if (defined $system_error);
            print_error(\%error);
        }
    }

    # return the split result
    return $split_result;

}


sub split_parser_options {

    my ($splitted_xml) = @_;

    # define index file name
    my $index_file_name_full_path = $splitted_xml->index_file_name();

    # open index file
    my $index_fh;
    unless (open($index_fh, '>', $index_file_name_full_path)) {
        my $system_error = $!;
	my ($filename, $dir, $ext) = fileparse($index_file_name_full_path, qr/\.[^.]*/);
	my %error;
    	$error{type}    = "ERROR";
	$error{message} = "cannot create index file";
	$error{message} .= " (" . $filename . $ext . ") ";
	$error{message} .= "in \"$dir\"";
	$error{message} .= ": $system_error" if (defined $system_error);
        print_error(\%error);
        return;
    }

    # initialize filehandles
    $splitted_xml->{index_fh}   = $index_fh;
    $splitted_xml->{current_fh} = $index_fh;

    # define handlers for start and end tags
    my $handlers = {Start   => sub { $splitted_xml->{handlers}->{Start}->($splitted_xml, @_); },
                    End     => sub { $splitted_xml->{handlers}->{End}->($splitted_xml, shift(@_)); },
                    Default => sub { $splitted_xml->{handlers}->{Default}->($splitted_xml, shift(@_)); },
                    XMLDecl => sub { split_parser_declaration->($splitted_xml, shift(@_)); }};
    
    return {Handlers => $handlers};

}


sub split_parser_start_tag {

    my ($splitted_xml, $split_parser, $element, %attributes) = @_;

    # remember current element and its attributes by pushing a hash onto the element stack
    push(@{$splitted_xml->{stack}}, { element_name=>$element, attributes=>{%attributes} });
   
    # prepare a new file when size is 0 (or when size is re-set) and the parser is deeper than the target depth
    if ( (! $splitted_xml->{current_size}) && ($split_parser->depth >= $splitted_xml->{target_depth})) {
    
        # increase sequence number
        $splitted_xml->{seq_number}++;

	# define file name
        my $file_name_full_path = $splitted_xml->file_name; 
	
        # open new splitted file
	my $out;
	unless (open($out, '>', $file_name_full_path)) {
            my $system_error = $!;
	    my ($filename, $dir, $ext) = fileparse($file_name_full_path, qr/\.[^.]*/);
	    my %error;
    	    $error{type}    = "ERROR";
	    $error{message} = "cannot create output file";
	    $error{message} .= " (" . $filename . $ext . ") ";
	    $error{message} .= "in \"$dir\"";
	    $error{message} .= ": $system_error" if (defined $system_error);
            print_error(\%error);
            die;
        }
        $splitted_xml->{current_fh} = $out if (defined $out);

        # make sure the XML declaration is defined
        $splitted_xml->{xml_declaration} = '' unless (defined $splitted_xml->{xml_declaration});

	# write XML declaration and header (if any) to new splitted file
	if ($splitted_xml->{current_fh}) {

            # print XML declaration to new splitted file
            print {$splitted_xml->{current_fh}} qq/$splitted_xml->{xml_declaration}\n\n/;

	    # define header filename (if any) and open it read-only 
	    my $in;
            if ((defined $splitted_xml->{header_dir}) && (defined $splitted_xml->{header_name}) && (defined $splitted_xml->{header_ext})) {
	        my $header_file_name_full_path = $splitted_xml->{header_dir} . $splitted_xml->{header_name} . $splitted_xml->{header_ext};
		unless (open($in, '<', $header_file_name_full_path)) {
                    my $system_error = $!;
	            my %error;
    	            $error{type}    = "WARNING";
	            $error{message} = "cannot open header file";
	            $error{message} .= " (" . $splitted_xml->{header_name} . $splitted_xml->{header_ext} . ") for reading";
	            $error{message} .= " in \"$splitted_xml->{header_dir}\"";
	            $error{message} .= ": $system_error" if (defined $system_error);
                    print_error(\%error);
                }
	    }

	    # print header to splitted file (skipping empty lines and XML declaration)
	    if (defined $in) {
	        
		# print new root start tag
		if (defined $splitted_xml->{new_root}) {
	            my $new_root_start_tag = '<' . $splitted_xml->{new_root} . '>';
                    print {$splitted_xml->{current_fh}} qq/$new_root_start_tag\n/;
		}

                # print header
	        while (<$in>) {
		    next if /^$/;
		    next if /^<\?xml(([\s\t]+version(.*)="([0-9\.]+)")?|([\s\t]+encoding="([a-zA-Z0-9-_]+)")?|([\s\t]+standalone="(yes|no)")?)*[\s\t]*\?>/;
		    print {$splitted_xml->{current_fh}} $_;
		}

		# close header
                unless (close($in)) {
                    my $system_error = $!;
	            my %error;
    	            $error{type}    = "WARNING";
	            $error{message} = "cannot close header file";
		    $error{message} .= " (" . $splitted_xml->{header_name} . $splitted_xml->{header_ext} . ")";
	            $error{message} .= " in \"$splitted_xml->{header_dir}\"";
	            $error{message} .= ": $system_error" if (defined $system_error);
                    print_error(\%error);
                }

	    }

	    # write current open element names and their attributes (if any)
	    my @context = $split_parser->context;
	    foreach my $element_name (@context) {
	        		
		# intialize start tag
		my $start_tag = '<' . $element_name;

		# check whether element has attributes and if so add them to the start tag
	        foreach my $element_from_stack (@{$splitted_xml->{stack}}) {
	            if (($element_name eq $element_from_stack->{element_name}) && ($element_from_stack->{attributes})) {
                        while ( my ($attribute, $attribute_value) = each( %{$element_from_stack->{attributes}} )) {
			    $start_tag .= ' ' . $attribute . '="' . $attribute_value . '"';
                        }
			last;
		    }
		}

		# finalize start tag
	        $start_tag .= '>';

		# write start tag to file
                if ( $start_tag !~ /Earth_Explorer_File/ || $splitted_xml->{data_file_ext} ne '.EEF') {
		    print {$splitted_xml->{current_fh}} qq/$start_tag\n/;
                }
	    }

	}

        # output splitted filename to index file
	if ($splitted_xml->{index_fh}) {
            print {$splitted_xml->{index_fh}} $splitted_xml->include($file_name_full_path);
	}

        $splitted_xml->{store_size} = 1;

    }

    $splitted_xml->{current_size} += length($split_parser->original_string) if ($splitted_xml->{store_size});
    print {$splitted_xml->{current_fh}} $split_parser->original_string if ($splitted_xml->{current_fh});    

}
  

sub split_parser_end_tag {
    
    my ($splitted_xml, $split_parser) = @_;

    # clean the element stack
    pop(@{$splitted_xml->{stack}});

    # calculate size
    $splitted_xml->{current_size} += length($split_parser->original_string) if ($splitted_xml->{store_size});

    # print to and, if needed, end file
    if (($splitted_xml->{current_size} <= $splitted_xml->{target_size}) 
	&& ($split_parser->depth < $splitted_xml->{target_depth})) {
	# end file (and change current filehandle to index file),
	# because depth is lower then target depth
	# (note: order [before print statement] is important)
        my @context = $split_parser->context;
	push @context, $split_parser->original_string;
        end_file_and_reset_size($splitted_xml, \@context);
    }
    if ( ($splitted_xml->{current_fh}) && (($splitted_xml->{current_fh} != $splitted_xml->{index_fh}) 
					   || ($split_parser->depth <= $splitted_xml->{target_depth}))) { 

        print {$splitted_xml->{current_fh}} $split_parser->original_string;

    }
    if (($splitted_xml->{current_size} > $splitted_xml->{target_size}) 
	&& ($split_parser->depth == $splitted_xml->{target_depth})) {
	# end file (and change current filehandle to index file),
	# because size is larger than target size
	# (note: order [after print statement] is important)
        my @context = $split_parser->context;
        end_file_and_reset_size($splitted_xml, \@context);
    }

}


sub end_file_and_reset_size {

    my ($splitted_xml, $context) = @_;
    
    if ($splitted_xml->{current_fh} != $splitted_xml->{index_fh}) {

        # close currently open elements (in reverse order)
        for (my $i = $#{$context}; $i >= 0; $i--) {
            my $element_name = $context->[$i];
	    my $end_tag      = ($element_name =~ /^<\/.*>$/) ? $element_name : '</' . $element_name . '>';
            if ($element_name !~ /Earth_Explorer_File/ || $splitted_xml->{data_file_ext} ne '.EEF') {
            print {$splitted_xml->{current_fh}} qq/$end_tag\n/ if (($splitted_xml->{current_fh}) && (defined $end_tag));
            }
        }

	# in case a header was added, print new root end tag
	if ((defined $splitted_xml->{header_dir}) && (defined $splitted_xml->{header_name}) && (defined $splitted_xml->{header_ext})) {
	    my $new_root_end_tag = '</' . $splitted_xml->{new_root} . '>' if (defined $splitted_xml->{new_root});
            print {$splitted_xml->{current_fh}} qq/$new_root_end_tag\n/ if (($splitted_xml->{current_fh}) && (defined $new_root_end_tag));
	}

	# close file
	if (defined $splitted_xml->{current_fh}) {
            unless (close($splitted_xml->{current_fh})) {
                my $system_error = $!;
	        my %error;
    	        $error{type}    = "WARNING";
	        $error{message} = "cannot close current open splitted file";
	        $error{message} .= ": $system_error" if (defined $system_error);
                print_error(\%error);
            }
        }

        # reset size
        $splitted_xml->{current_size} = 0; 
        $splitted_xml->{store_size}   = 0;
	
        $splitted_xml->{current_fh} = $splitted_xml->{index_fh};

    }

}


sub split_parser_default {

    my ($splitted_xml, $split_parser) = @_;

    if ($splitted_xml->{store_size}) {
        
	$splitted_xml->{current_size} += length($split_parser->original_string);
        
	if ($split_parser->depth < $splitted_xml->{target_depth}) {
	    my @context = $split_parser->context;
	    end_file_and_reset_size($splitted_xml, \@context);
	}

    }

    print {$splitted_xml->{current_fh}} $split_parser->original_string if ($splitted_xml->{current_fh});

}


sub split_parser_declaration {

    my ($splitted_xml, $split_parser) = @_;

    $splitted_xml->{xml_declaration} = $split_parser->original_string || '';

    if (defined $splitted_xml->{index_fh}) {
        print {$splitted_xml->{index_fh}} qq/$splitted_xml->{xml_declaration}\n\n/ unless ($splitted_xml->{seq_number});
    }
       
    if (( ! $splitted_xml->{xml_declaration}) || ($splitted_xml->{xml_declaration} =~ /encoding\s*=\s*["']utf-?8["']/i)) {
        $splitted_xml->{utf8_encoded} = 1;
        $split_parser->setHandlers(Char => \&split_parser_contents);
    }

}


sub split_parser_contents {

    my ($splitted_xml, $split_parser) = (shift, shift);

    print {$splitted_xml->{current_fh}} $_[0] if ($splitted_xml->{current_fh});

}


sub parse_index_file {

    my $index_file_full_path = shift;

    # declare return variable
    my $splitted_filenames;

    # return if we have nothing to parse
    return unless (defined $index_file_full_path);

    # define a new XML parser object with options
    my $index_parser_options = index_parser_options();
    my $index_parser         = XML::Parser->new(%{$index_parser_options}) if (defined $index_parser_options);

    # parse and split the data file
    $splitted_filenames = $index_parser->parsefile($index_file_full_path);

    # return splitted filenames
    return $splitted_filenames;

}


sub index_parser_options {

    # define handlers
    my $handlers = {Start => \&index_parser_start_tag,
                    Final => \&index_parser_results};
    
    return {Handlers => $handlers};

}


sub index_parser_start_tag {

    my ($index_parser, $element, %attributes) = @_;

    if (($element eq 'Split') && (exists $attributes{filename})) {
        
	# store the filename in a variable
	my $filename_full_path = $attributes{filename};
	
	# define the number of currently open tags 
	my @open_elements = $index_parser->context;
	my $number_of_open_elements = $#open_elements + 1;
	
	# define main product and, if available, the sub product
	my $main_product = $open_elements[1];
	my $sub_product  = $open_elements[2] if ((defined $main_product) && ($number_of_open_elements == 4));

	#For some L1b products, we do not extract all the data
        if ( ($filename_full_path =~ /EGG_NOM_1/ && $main_product !~ /EGG_CCD|EGG_GGT|EGG_IAQ/) || 
	     ($filename_full_path =~ /SST_NOM_1/ && $main_product !~ /SST_COV|SST_PVT/)  || 
	     ($filename_full_path =~ /STR_VC/    && $main_product !~ /STR_VC/)   ) {	   
	   return
	}

        # store filename in index_parser, sorted by (sub)product
        if (defined $filename_full_path) {
	    if (defined $sub_product) {
	        push @{$index_parser->{splitted_filenames}->{$sub_product}}, $filename_full_path;
	    } elsif (defined $main_product) {
	        push @{$index_parser->{splitted_filenames}->{$main_product}}, $filename_full_path;
	    }
	}

    } 

}


sub index_parser_results {

    my ($index_parser) = @_;

    # declare return variable
    my %splitted_filenames;

    # store the splitted filenames sort by their (sub)product in a hash
    if (exists $index_parser->{splitted_filenames}) {
        %splitted_filenames = %{$index_parser->{splitted_filenames}};
    }

    # return a reference to the (sorted) splitted filenames
    return \%splitted_filenames;

}


sub transform_files {

    my ($arguments) = @_;

    my $filenames                = $arguments->{filenames};
    my $input_product_instance   = $arguments->{input_product_inst};
    my $stylesheet_doc           = $arguments->{stylesheet};
    my $mode_ids                 = $arguments->{mode_ids};
    my $max_recursion_depth      = $arguments->{max_recusrion_depth};
    my $output_ext               = $arguments->{output_ext};
    my $output_header_dir        = $arguments->{output_header_dir};
    my $output_data_file_dir     = $arguments->{output_data_file_dir};
    my $mission                  = $arguments->{mission};
    my $file_class               = $arguments->{file_class};
    my $validity_period          = $arguments->{validity_period};
    my $file_version             = $arguments->{file_version};

    # define mode ID's
    return unless ((defined $mode_ids) && (ref $mode_ids eq "HASH"));
    my $header                   = $mode_ids->{'header'};
    my $header_and_last_file     = $mode_ids->{'header_last_file'};
    my $data_block               = $mode_ids->{'data_block'};
    my $data_block_and_last_file = $mode_ids->{'data_block_last_file'};
    return unless ((defined $header) && (defined $header_and_last_file) && (defined $data_block) && (defined $data_block_and_last_file));

    # make sure we have stylesheets
    return unless ((defined $stylesheet_doc) && (ref $stylesheet_doc eq "HASH"));
    return unless ((exists $stylesheet_doc->{$header}) && (exists $stylesheet_doc->{$data_block}));
    return unless ((defined $stylesheet_doc->{$header}) && (defined $stylesheet_doc->{$data_block}));

    # define new XSLT object for transforming files
    # and set its options (maximum recursion depth)
    my $xslt = XML::LibXSLT->new();
    $xslt->max_depth($max_recursion_depth);

    # define a new parser object for parsing stylesheets
    my $stylesheet_parser = XML::LibXML->new();

    # parse style sheets
    my (%parsed_stylesheet_doc, %parsed_stylesheet);
    for my $stylesheet_type ($header, $data_block) {
	unless (file_is_not_ok($stylesheet_doc->{$stylesheet_type}, 'xsl')) {
            $parsed_stylesheet_doc{$stylesheet_type} = $stylesheet_parser->parse_file($stylesheet_doc->{$stylesheet_type});
            if (defined $parsed_stylesheet_doc{$stylesheet_type}) {
                eval { $parsed_stylesheet{$stylesheet_type} = $xslt->parse_stylesheet($parsed_stylesheet_doc{$stylesheet_type}); };
                if ($@) {
                    my $stylesheet_error = $@;
	            if (defined $stylesheet_error) {
	                my $stylesheet_error_max_line_length = 160;
		        chomp $stylesheet_error;
		        $stylesheet_error =~ s/\.*\s*\n/ /g if (length($stylesheet_error) <= $stylesheet_error_max_line_length);
	            }
	            my %error;
	            $error{type}     = "ERROR";
	            $error{message}  = "in stylesheet";
	            $error{message} .= " ($stylesheet_error)" if ((defined $stylesheet_error) && ($stylesheet_error ne ''));
	            $error{message} .= ", $stylesheet_type will not be transformed";
	            print_error(\%error);
	        }
            }
	}

    }
    
    # define header output file and open it
    my $output_header_fh;
    if ((exists $parsed_stylesheet{$header}) && (defined $parsed_stylesheet{$header})) {
    
        # define output header file
	my $output_header_full_path;
        if (exists $output_ext->{+HEADER}) {
            $output_header_full_path = $output_header_dir . $g_dir_separator . $input_product_instance . $output_ext->{+HEADER};
	    $output_header_full_path =~ s/$g_multiple_dir_separators/$g_dir_separator/g;
        }

        # if desired, print some extra information before transforming
        print STDOUT "OUTPUT HEADER   : $output_header_full_path\n" if ($g_verbose);

        # open output header file
	if (defined $output_header_full_path) {
	    unless (open($output_header_fh, '>', $output_header_full_path)) {
	        my $system_error = $!;
		my ($filename, $dir, $ext) = fileparse($output_header_full_path, qr/\.[^.]*/);
		my %error;
		$error{type}    = "ERROR";
		$error{message} = "cannot create output header file";
		$error{message} .= " (" . $filename . $ext . ") ";
		$error{message} .= "in \"$dir\"";
		$error{message} .= ": $system_error" if (defined $system_error);
		print_error(\%error);
	    }
	}

    }

    # transform header and/or data
    if (((exists $parsed_stylesheet{$data_block}) && (defined $parsed_stylesheet{$data_block})) ||
        ((exists $parsed_stylesheet{$header}) && (defined $parsed_stylesheet{$header}))) {

        # initialize transform warning flags
        my $header_transform_warning = 0;
        my $data_transform_warning   = 0;

        # set first time flag (for next loop)
	my $first_output_product = 1;

        OUTPUT_PRODUCT:
        foreach my $output_product (keys %{$filenames}) {

            # initialize (or reset) Base64
	    # encoded and decoded flags
	    my $base64_encoded = 0;
	    my $base64_decoded = 0;

	    # declare (or reset) output file name and handles
	    my ($output_data_file_full_path, $output_data_fh, $base64_decoded_fh);

            if ((exists $parsed_stylesheet{$data_block}) && (defined $parsed_stylesheet{$data_block})) {

	        # define output data file
                if (exists $output_ext->{+DATA}) {
                    my $product_name_length        = 10;
                    my $padding_character          = '_';
                    my $formatted_output_product   = $output_product . ($padding_character x ($product_name_length - length($output_product)));
	            my $output_product_instance    = $mission . '_' . $file_class . '_' . $formatted_output_product . '_' . $validity_period . '_' . $file_version;
		    $base64_encoded                = 'base64_encoded' if (substr($output_product, 5, 2) =~ /^RP$/io);
		    my $data_file_output_extension = ($base64_encoded) ? substr($output_ext->{+DATA}, 0, 1) . $base64_encoded . '_'. substr($output_ext->{+DATA}, 1) : $output_ext->{+DATA};
                    $output_data_file_full_path    = $output_data_file_dir . $g_dir_separator . $output_product_instance . $data_file_output_extension;
	            $output_data_file_full_path    =~ s/$g_multiple_dir_separators/$g_dir_separator/g;
                }

	        # open output data file and, if desired, print some extra information
	        if (defined $output_data_file_full_path) {

                    if ($g_verbose) {
		    	my $stdout_output_data_file_full_path = $output_data_file_full_path;
		        ($stdout_output_data_file_full_path = $output_data_file_full_path) =~ s/\.${base64_encoded}_/\./ if ($base64_encoded);
                        print STDOUT "OUTPUT DATA FILE: $stdout_output_data_file_full_path\n";
                    }

		    my $open_mode = ($base64_encoded) ? '+>' : '>'; 
	            unless (open($output_data_fh, $open_mode, $output_data_file_full_path)) {
	                my $system_error = $!;
		        my ($filename, $dir, $ext) = fileparse($output_data_file_full_path, qr/\.[^.]*/);
		        my %error;
		        $error{type}    = "ERROR";
		        $error{message} = "cannot create output data file";
		        $error{message} .= " (" . $filename . $ext . ") ";
		        $error{message} .= "in \"$dir\"";
		        $error{message} .= ": $system_error" if (defined $system_error);
		        print_error(\%error);
		        last OUTPUT_PRODUCT;
	            }

                    print STDOUT "Transforming splitted input data files. This can take a while, please wait...\r" if ($g_verbose);

	        } else {

	            last OUTPUT_PRODUCT;

	        }

            }

            # set first and last filename flag (for next loop)
	    my $first_filename = 1;
	    my $last_filename  = 0;

	    # determine the last input filename for the current output product, needed for setting the last filename flag
	    my $last_input_filename_full_path = $filenames->{$output_product}->[$#{$filenames->{$output_product}}];

            # parse files and transform
	    FILENAME:
	    foreach my $input_filename_full_path (@{$filenames->{$output_product}}) {

                # set last time flag
		$last_filename = 1 if ($input_filename_full_path eq $last_input_filename_full_path);

	        # define the XSLT mode, parsed to the stylesheet as parameter,
	        # to avoid printing the data header more than once and to determine
		# whether it's the last input file (for writing end of file markers)
	        my $mode;
		if ($first_filename && $last_filename) {
		    $mode = $header_and_last_file;
		} elsif ($first_filename) {
		    $mode = $header
		} elsif ($last_filename) {
		    $mode = $data_block_and_last_file;
		} else {
		    $mode = $data_block;
		}

                # define new LibXML object and parse the data
                my $parser      = XML::LibXML->new();
	        my $parsed_data = $parser->parse_file($input_filename_full_path);
		unless (defined $parsed_data) {
		    $first_filename = 0;
		    last FILENAME;
		}

                # transform header and output to filehandle; only once and close output header file immediately afterwards
                if (($first_output_product) && ($first_filename) && (exists $parsed_stylesheet{$header}) && (defined $parsed_stylesheet{$header})) {
                    my $transform_header_results = transform_parsed_data($parsed_data, $parsed_stylesheet{$header});
                    if ((defined $output_header_fh) && (defined $transform_header_results)) {
		        $parsed_stylesheet{$header}->output_fh($transform_header_results, $output_header_fh);
		    } else {
                        $header_transform_warning = 1;
	            }
		    if (defined $output_header_fh) {
                        unless (close($output_header_fh)) {
                            my $system_error = $!;
	                    my %error;
    	                    $error{type}    = "WARNING";
	                    $error{message} = "cannot close output header file";
	                    $error{message} .= ": $system_error" if (defined $system_error);
                            print_error(\%error);
                        }
                    }
                    last OUTPUT_PRODUCT unless ((exists $parsed_stylesheet{$data_block}) && (defined $parsed_stylesheet{$data_block}));
		}

                # transform data and output to filehandle
		if ((exists $parsed_stylesheet{$data_block}) && (defined $parsed_stylesheet{$data_block})) {

                    # transform data
		    my $transform_data_results = transform_parsed_data($parsed_data, $parsed_stylesheet{$data_block}, $mode, $output_product);

                    if ((defined $output_data_fh) && (defined $transform_data_results)) {

		        # write to output file
			$parsed_stylesheet{$data_block}->output_fh($transform_data_results, $output_data_fh);
			
			if ($base64_encoded) {

			    # define the file to decode to
			    my $base64_decoded_file_full_path;
			    ($base64_decoded_file_full_path = $output_data_file_full_path) =~ s/\.${base64_encoded}_/\./;

			    # define filehandles and open output file
			    my $base64_encoded_fh = $output_data_fh;
			    unless (open($base64_decoded_fh, '>', $base64_decoded_file_full_path)) {
	                        my $system_error = $!;
		                my ($filename, $dir, $ext) = fileparse($base64_decoded_file_full_path, qr/\.[^.]*/);
		                my %error;
		                $error{type}     = "WARNING";
		                $error{message}  = "cannot create output data file";
		                $error{message} .= " (" . $filename . $ext . ")";
		                $error{message} .= " in \"$dir\"";
		                $error{message} .= " for writing decoded Base64 encoded data";
		                $error{message} .= ": $system_error" if (defined $system_error);
		                print_error(\%error);
	                    }

			    # decode Base64 encoded data
			    if ((defined $base64_encoded_fh) && (defined $base64_decoded_fh)) {
			        $base64_decoded = decode_base64_data({input_fh  => $base64_encoded_fh,
			                                              output_fh => $base64_decoded_fh});
			    }

			}

		    }

		}

                # cleanup of temporary splitted data file
		if (($g_cleanup) && (-e $input_filename_full_path)) {
		    unless ($header_transform_warning || $data_transform_warning) {
		        unless (unlink $input_filename_full_path) {
	                    my $system_error = $!;
		            my ($filename, $dir, $ext) = fileparse($input_filename_full_path, qr/\.[^.]*/);
		            my %error;
		            $error{type}    = "WARNING";
		            $error{message} = "cannot delete (temporary) splitted XML data file";
		            $error{message} .= " (" . $filename . $ext . ")";
		            $error{message} .= " from \"$dir\"";
		            $error{message} .= ": $system_error" if (defined $system_error);
		            print_error(\%error);
			}
                    }
		}

		# unset first and last time flags
	        $first_filename = 0;
		$last_filename  = 0;

            }

            # close output data file(s)
            if (defined $output_data_fh) {
                unless (close($output_data_fh)) {
                    my $system_error = $!;
	            my %error;
    	            $error{type}    = "WARNING";
	            $error{message} = "cannot close output data file";
	            $error{message} .= ": $system_error" if (defined $system_error);
                    print_error(\%error);
                }
            }
	    if (defined $base64_decoded_fh) {
	        unless (close($base64_decoded_fh)) {
                    my $system_error = $!;
	            my %error;
    	            $error{type}    = "WARNING";
	            $error{message} = "cannot close Base64 decoded output data file";
	            $error{message} .= ": $system_error" if (defined $system_error);
                    print_error(\%error);
                }
		if ($base64_decoded && $g_cleanup) {
		    unless (unlink $output_data_file_full_path) {
	                my $system_error = $!;
		        my ($filename, $dir, $ext) = fileparse($output_data_file_full_path, qr/\.[^.]*/);
		        my %error;
		        $error{type}    = "WARNING";
		        $error{message} = "cannot delete Base64 encoded input data file";
		        $error{message} .= " (" . $filename . $ext . ")";
		        $error{message} .= " from \"$dir\"";
		        $error{message} .= " after successful decoding";
		        $error{message} .= ": $system_error" if (defined $system_error);
		        print_error(\%error);
		    }
		}

	    }

	    # unset first time flag
	    $first_output_product = 0;

	}

        # if desired, print some extra information after transforming
	# (extra spaces at the end are needed to clear previous message)
	print STDOUT "Finished transforming splitted input data files.                             \n" if ($g_verbose);

	# return OK (=true) unless there was a warning
	return 'OK' unless ($header_transform_warning || $data_transform_warning);

    }

    return;

}


sub transform_parsed_data {

    my $data       = shift;
    my $stylesheet = shift;
    my $mode       = shift;
    my $product    = shift;

    # make sure we have something to transform
    return unless ((defined $data) && (defined $stylesheet));

    # define return variable
    my $transform_results;

    # transform
    if ((defined $mode) && (defined $product)) {
        eval { $transform_results = $stylesheet->transform($data, XML::LibXSLT::xpath_to_string(Mode => "$mode"), XML::LibXSLT::xpath_to_string(Product => "$product")); };
    } elsif (defined $mode) {
        eval { $transform_results = $stylesheet->transform($data, XML::LibXSLT::xpath_to_string(Mode => "$mode")); };
    } elsif (defined $product) {
        eval { $transform_results = $stylesheet->transform($data, XML::LibXSLT::xpath_to_string(Product => "$product")); };
    } else {
        eval { $transform_results = $stylesheet->transform($data); };
    }
    
    # error handling and return
    if ($@) {
        
	my $transform_warning = $@;
	if (defined $transform_warning) {
	    my $transform_warning_max_line_length = 160;
	    chomp $transform_warning;
	    $transform_warning =~ s/\.*\s*\n/ /g if (length($transform_warning) <= $transform_warning_max_line_length);
	}
	
	my %error;
	$error{type}     = "WARNING";
	$error{message}  = "in transforming data";
	$error{message} .= " ($transform_warning)" if ((defined $transform_warning) && ($transform_warning ne ''));
	print_error(\%error);

	return;

    } else {

        return $transform_results;

    }

}


sub decode_base64_data {

    my ($arguments) = @_;
    
    my $input_fh  = $arguments->{'input_fh'};
    my $output_fh = $arguments->{'output_fh'};

    # rewind the input and output filehandels
    # (just to be sure to start at the beginning)
    seek $input_fh, 0, 0;
    seek $output_fh, 0, 0;
    
    # declare variable for read buffer
    my $read_buffer;
		    
    # put the output file in binary mode and decode
    binmode($output_fh);
    while ($read_buffer = readline($input_fh)) {
        print $output_fh decode_base64($read_buffer);
    }

    return 1;

}


sub print_error {

    my $error = shift;

    # set a minimal error message length to overwrite previous
    # log messages, if needed (determined by a \r character)
    my $minimal_error_message_length = 114;

    # print final error message to standard error
    if ((defined $error) && (ref $error eq "HASH")) {
        my $final_error_message = $error->{type} . ': ' . $error->{message};
	if (length($final_error_message) < $minimal_error_message_length) {
	    printf STDERR "%-${minimal_error_message_length}s\n", $final_error_message;
	} else {
            print STDERR "$error->{type}: $error->{message}\n";
	}
    } else {
	printf STDERR "%-${minimal_error_message_length}s\n", "ERROR: unknown error";
    }

    return;

}


sub check_command_line_filenames {

    my $filenames_from_cmd_line = shift;
    my $extension               = shift;
    my $products                = shift;

    # declare return variables
    my (%filenames, %header_found, %datafile_found);

    # go through the filenames from the command line
    CMD_LINE_FILE:
    foreach my $filename (@{$filenames_from_cmd_line}) {

        # initialize/reset duplicate filename flag
	my $duplicate_name_error = 0;

        # split in directory, filename and extension
	# (note that due to the fileparse function,
	# the return variables will always be defined)
        my ($product_instance, $dir, $ext) = fileparse($filename, qr/\.[^.]*/);
	
	# check that a Level 1b product has an EEF extension
	if (($product_instance =~ /(.*)_1b_(.*)/) && $ext ne '.EEF') {
	    my %error;
	    $error{type}    = "ERROR ";
	    $error{message} = " The L1b product $filename must have an EEF extension";
            print_error(\%error);
	    exit;
        }
	    

        # check whether the filename from the command
	# line contains one of the allowed products
	my $product_allowed = 0;
	foreach my $allowed_product (keys %{$products}) {
	    $product_allowed = 1 if ($product_instance =~ /$allowed_product/i);
	    last if ($product_allowed);
	}
	unless ($product_allowed) {
	    my %error;
	    $error{type}    = "WARNING";
	    $error{message} = "\"$filename\" will not be transformed (product not allowed)";
            print_error(\%error);
	    next CMD_LINE_FILE;
	}

        # store file info in filenames hash and remember
	# whether a header and a data file were found
	if (exists $extension->{+INPUT}) {
	    if ((exists $extension->{+INPUT}->{+HEADER}) && ($ext eq $extension->{+INPUT}->{+HEADER})) {
                $header_found{$product_instance} = 1;
		if (exists $filenames{$product_instance}{+HEADER}{ext}) {
	            $duplicate_name_error = 1;
		} else {
		    $filenames{$product_instance}{+HEADER}{ext} = $ext;
		    $filenames{$product_instance}{+HEADER}{dir} = $dir;
		}
	    } elsif ((exists $extension->{+INPUT}->{+DATA}) && ($ext eq $extension->{+INPUT}->{+DATA})) {
                $datafile_found{$product_instance} = 1;
		if (exists $filenames{$product_instance}{+DATA}{ext}) {
	            $duplicate_name_error = 1;
		} else {
		    $filenames{$product_instance}{+DATA}{ext} = $ext;
		    $filenames{$product_instance}{+DATA}{dir} = $dir;
		}
	    } elsif ($ext eq '.EEF') {
               $datafile_found{$product_instance} = 1;
               $filenames{$product_instance}{+DATA}{ext} = $ext;
               $filenames{$product_instance}{+DATA}{dir} = $dir;
               # overwrite the configuration 
               $extension->{+INPUT}->{+HEADER} = '';
               $extension->{+INPUT}->{+DATA}   = '.EEF'; 
            }  
	}

	if ($duplicate_name_error) {
	    my %error;
	    $error{type}    = "WARNING";
	    $error{message} = "filenames with the same name encountered on the command line ($filename), only the first one is used";
            print_error(\%error);
	}

    }

    return (\%filenames, \%header_found, \%datafile_found);

}


sub parse_eef_filename {

    my $filename = shift;
 
    return unless (defined $filename);

    # strip the path from the filename
    my $file_basename = basename($filename);

    # match the file basename according to the EEF naming convention
    $file_basename =~ /^(\w{2})_(\w{4})_(\w{10})_(\w{31})_(\d{4})\.(EEF|HDR|DBL|TAR|TGZ|ZIP|GZ)$/io;

    my ($mission, $file_class, $product, $productinstance, $file_version, $product_type) = ($1, $2, $3, $4, $5, $6);

    # if one of the variables is undefined return 'undefined' 
    unless ((defined $mission) &&
            (defined $file_class) &&
	    (defined $product) &&
	    (defined $productinstance) &&
	    (defined $file_version) &&
	    (defined $product_type)) {
	my %error;
	$error{type}    = "ERROR";
        $error{message} = "\"$file_basename\" does not follow the EEF naming convention";
        print_error(\%error);
	return;
    }

    # remove trailing underscores
    $product         =~ s/_+$//o;
    $productinstance =~ s/_+$//o;

    # make sure the mission, the file class, the product
    # and the product type are in upper case
    $mission      =~ tr/a-z/A-Z/;
    $file_class   =~ tr/a-z/A-Z/;
    $product      =~ tr/a-z/A-Z/ if $product !~ /1b/;
    $product_type =~ tr/a-z/A-Z/;
    
    # check two more EEF naming conventions for the GOCE satellite mission
    my $mission_full_name  = 'GOCE satellite mission';
    my $mission_identifier = substr($mission_full_name, 0, 2);
    if ( $mission !~ /^${mission_identifier}$/ ) {
        # log error and return 'undefined'
	my %error;
	$error{type}     = "ERROR";
        $error{message}  = "\"$file_basename\" does not belong to the $mission_full_name";
	$error{message} .= "\n                                all files should start with the mission identifier \"$mission_identifier\"";
        print_error(\%error);
        return;
    } elsif ( $file_class !~ /^(?:OPER|CONS|RPRO|TEST|TD\d{2})$/o) {
        # log error and return 'undefined'
	my %error;
	$error{type}     = "ERROR";
        $error{message}  = "\"$file_basename\" does not follow the EEF naming convention";
        $error{message} .= "\n                                file class \"$file_class\" is not allowed";
        print_error(\%error);
        return;
    }

    # verify the product name and
    # determine the level of the product
    # based on its product name
    my $level;
    if ($product =~ /^\w{3}_\w{3}_((?:1|2)(?:b|B|C|I)?)$/o) {
        $level = $1;
    } elsif ($product =~ /^(?:AUX|ANC|MPL)_\w{3,6}$/o) {
	$level = 'AUX';
	if ($product =~ /_(1|2)(?:b|B|C|I)?$/o) {
            $level .= $1;
	} else {
	    $level .= '1';
	}
    } else {
        # log error and return 'undefined'
	my %error;
	$error{type}     = "ERROR";
        $error{message}  = "\"$file_basename\" does not follow the EEF naming convention";
        $error{message} .= "\n                                unknown product type ($product)";
        print_error(\%error);
        return;
    }

    return ($mission, $file_class, $product, $level, $productinstance, $file_version, $product_type);

}


sub trim {

    my @out = @_;

    for (@out) {

        # remove leading and trailing whitespace
	s/^\s+//o; s/\s+$//o;

    }

    return wantarray ? @out : $out[0];

}


sub dir_is_not_ok {

    my $dir  = shift;
    my $type = shift;

    # make sure we have something to check
    return 1 unless (defined $dir);

    # initialize error flag
    my $error_encountered = 0;

    # initialize error type and message
    my %error;
    $error{type}     = "ERROR";
    $error{message}  = "$type " if (defined $type);
    $error{message} .= "directory \"$dir\" ";

    # check whether directory exists
    unless (-e $dir) {
        $error{message} .= "does not exist";
        $error_encountered = 1;
    }
    
    # check whether it is a directory
    if ((! $error_encountered) && (! -d $dir)) {
        $error{message} .= "is not a directory";
        $error_encountered = 1;
    }

    # check whether it is writable
    if ((! $error_encountered) && (! -w $dir )) {
        $error{message} .= "is not writable" unless ($error_encountered);
        $error_encountered = 1;
    }

    # return 0 or 1
    if ($error_encountered) {
        print_error(\%error);
        return 1;
    } else {
        return 0;
    }

}


sub file_is_not_ok {

    my $file = shift;
    my $type = shift;

    # make sure we have something to check
    return 1 unless (defined $file);

    # initialize error flag
    my $error_encountered = 0;

    # initialize error type and message
    my %error;
    $error{type}     = "ERROR";
    $error{message}  = "$type " if (defined $type);
    $error{message} .= "file \"$file\" ";

    # check whether file exists
    unless ( -e $file ) {
        $error{message} .= "does not exist";
        $error_encountered = 1;
    }

    # check whether it is readable
    if ((! $error_encountered) && (! -r $file)) {
        $error{message} .= "is not readable";
        $error_encountered = 1;
    }

    # check whether it is a text file (not binary)
    if ((! $error_encountered) && (! -T $file )) {
        $error{message} .= "is not a text file";
        $error_encountered = 1;
    }

    # return 0 or 1
    if ($error_encountered) {
        print_error(\%error);
        return 1;
    } else {
        return 0;
    }

}


sub usage {

    my $default_config_file = shift;

    print STDERR <<EndOfUsage;
Usage:
   $g_program [--help] [--verbose] [--nocleanup] [--config <config file>] [--ouput <output dir>] <EEF file(s)>
Where:
   --help			show usage and version information and exit (all other options are ignored)
   --verbose                    show extra information during processing
   --nocleanup                  do not cleanup temporary splitted files and their index file   
   --config			specify a different XML configuration file (default: $default_config_file)
   --output <output dir>	specify a different output directory (default: same as input files)
   <EEF file(s)>		is one or more input files (separated by spaces) in EEF format which need to be transformed
Version:
   $g_version

EndOfUsage
    exit 1;
    
}




package xml_split::splitted_xml;

sub new {

    my ($ref, $options) = @_;

    my $splitted_xml = bless $options, $ref;
    $splitted_xml->{seq_number} = 0;

    return $splitted_xml;

}


sub file_name {

    my ($splitted_xml) = @_;

    my $number    = sprintf("%0$splitted_xml->{digits}d", $splitted_xml->{seq_number});
    my $file_name = "$splitted_xml->{output_dir}$splitted_xml->{data_file_name}-$number$splitted_xml->{ext}";

    return $file_name;

}


sub index_file_name {

    my ($splitted_xml) = @_;

    my $number    = sprintf("%0$splitted_xml->{digits}d", 0);
    my $file_name = "$splitted_xml->{output_dir}$splitted_xml->{data_file_name}-$number$splitted_xml->{ext}";

    return $file_name;

}

1;




package xml_split::splitted_xml::split_parser;

import xml_split::splitted_xml;
use base 'xml_split::splitted_xml';

sub include {

    my ($splitted_xml, $file_name) = @_;
    
    return qq/<Split filename="$file_name"\/>/;

}

1;
