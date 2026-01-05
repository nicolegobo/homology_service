
=head1 NAME

HomologySearch - Homology search application.

=head1 DESCRIPTION

The homology search application searches one or more databases for
occurrences of sequences in its input set. The application parameters
contain all the values we describe below.

=head2 Database types

The search databases contain sequences of four different types. Each of these
types is given a unique suffix in the name for the databases containing
data of the given type. This suffix is also used as the key to select a 
database type in the C<db_type> parameter.

=over 4

=item * 

Protein features. These are the amino acid sequences of the protein encoding genes in 
the genome. The suffix is I<faa>.

=item *

DNA features. These are the DNA-based sequences of any genes as called in the genome. 
The suffix is I<ffn>.

=item * 

RNA features. The suffix is I<frn>.

=item *

Genome sequences (contigs). The suffix is I<fna>.

=back

=head2 Specifying input sequences

The type of input (DNA or amino acid) is specified by the C<input_type> parameter. 
It has one of the following values:

=over 4

=item dna

Input is DNA sequences.

=item aa

Input is amino acid sequences.

=back

The source of input is specified by the C<input_source> parameter. It has one of the
following values:

=over 4

=item id_list

Use a set of sequences from the PATRIC database named by the parameter C<input_id_list>. 

=item fasta_data

Use the sequence data from the parameter C<input_fasta_data>. 

=item fasta_file

Use the sequence data from the workspace file at path C<input_fasta_file>.

=item feature_group

Use a set of sequences from the PATRIC database named by the identifiers in the feature group
whose path is defined by the parameter C<input_feature_group>.
    
=back

=head2 Specifying search database

The source of the database to search is specified by the C<db_source> parameter. It
has one of the following values:

=over 4

=item fasta_data

Use the sequence data from the parameter C<db_fasta_data>.

=item fasta_file

Use the sequence data from the workspace file at path C<db_fasta_file>.

=item genome_list

Use the data from the set of PATRIC genome ids specified in C<db_genome_list>.

=item taxon_list

Use the data from genomes in PATRIC at or below the given taxon identifiers. given in C<db_taxon_list>.

=item precomputed_database

Use the precomputed database specified by C<db_precomputed_database>. The list of available databases
may be queried using the HomologySearch service.

=back

=head2 Specifying behavior

The BLAST program to use is defined by the C<blast_program> parameter. Valid values
are C<blastp>, C<blastn>, C<blastx>, C<tblastn>, C<tblastx>. A value is not required; if not
specified, the appropriate program is chosen for the given input and database types.
For DNA to DNA searches, C<blastn> is the default.

=head1 PRECOMPUTED DATABASES

To accelerate common searches, we maintain a number of precomputed BLAST databases:

=over 4

=item * 

Each of the 26 large and/or pathogen specific bacterial genera has a BLAST database containing all genomes in the genus.

=item * 

Each of the 21 viral families has a BLAST database containing all genomes in the family.

=item * 

A database of all the reference and representative genomes in the bacterial and archaeal taxonomies (~6700 genomes).

=back




=head1 METHODS

=cut

package Bio::P3::HomologySearch::HomologySearch;

use Bio::P3::HomologySearch::Config qw(blast_db_search_path blast_sqlite_db);
use Bio::P3::HomologySearch::BlastDatabases;
use Bio::P3::HomologySearch::BlastDatabasesSQL;

use P3DataAPI;
use gjoseqlib;
use strict;
use Data::Dumper;
use POSIX;
use Cwd;
use base 'Class::Accessor';
use gjoseqlib;
use JSON::XS;
use Module::Metadata;
use List::Util qw(any none);
use IPC::Run qw(run start finish);
use File::Temp qw(:seekable);
use File::Basename;
use File::Copy qw(copy);
use File::Path qw(make_path remove_tree);
use Statistics::Descriptive;
use Template;
use Text::CSV_XS qw(csv);
use Bio::KBase::AppService::ClientExt;
use Bio::KBase::AppService::AppConfig qw(data_api_url);
use Bio::KBase::AppService::FastaParser 'parse_fasta';
use Bio::P3::Workspace::WorkspaceClientExt;

use File::Slurp;

#
# Configure PATH to point at the correct location for BLAST binaries.
#
# If we have $KB_TOP/services/homology_service/bin/blastp, that is our
# deployed path and use that.
#
# Otherwise, assume we are in a dev container where we will look in
# $KB_TOP/modules/homology_service/blast.bin.
#


{
    my @blast_paths = ("$ENV{KB_TOP}/services/homology_service/bin",
		       "$ENV{KB_TOP}/modules/homology_service/blast.bin");
    my @found = grep { -x "$_/blastp" } @blast_paths;
    if (@found)
    {
	print STDERR "Initializing BLAST path to $found[0]\n";
	$ENV{PATH} = "$found[0]:$ENV{PATH}";
    }
    else
    {
	print STDERR "No internal BLAST found, using default path $ENV{PATH}\n";
    }
    my $out;
    run(["blastp", "-version-full"], ">", \$out);
    print STDERR "BLAST version: $out\n";
}

__PACKAGE__->mk_accessors(qw(app app_def params token task_id
			     blast_program blast_params
			     blast_sqlite_db
			     work_dir stage_dir
			     output_base output_folder 
			     contigs app_params  api
			     database_path
			     blastdbs
			     json
			    ));


#
# $search_type_map{input-type}->{db-type} = valid-blast-programs
#
our %search_type_map =
    (aa => {
	faa => ['blastp'],
	ffn => ['tblastn'],
	frn => ['tblastn'],
	fna => ['tblastn'],
    },
     dna => {
	 faa => ['blastx'],
	 ffn => ['blastn', 'tblastx'],
	 frn => ['blastn', 'tblastx'],
	 fna => ['blastn', 'tblastx'],
     });
				    
our %makeblastdb_dbtype = (faa => "prot",
			   ffn => "nucl",
			   frn => "nucl",
			   fna => "nucl");

sub new
{
    my($class) = @_;

    my $api = P3DataAPI->new(data_api_url);
    my $self = {
	api => $api,
	database_path => blast_db_search_path,
	json => JSON::XS->new->pretty,
	short_feature_threshold => 30,
    };

    $self->{blast_sqlite_db} = blast_sqlite_db;

    if ($self->{blast_sqlite_db})
    {
	$self->{blastdbs} = Bio::P3::HomologySearch::BlastDatabasesSQL->new(blast_db_search_path,
									    $self->{blast_sqlite_db}, $api);
    }
    else
    {
	$self->{blastdbs} = Bio::P3::HomologySearch::BlastDatabases->new(blast_db_search_path, $api);
    }

    bless $self, $class;

    return $self;
}

=head2 preflight

Compute CPU & estimated time.

We make really gross general estimates on the runtime here by computing the size in bases
of the input (we can do this fairly accurately since we're pulling input from immediate data
or from files which we can do a C<stat()> on to get size from.

We also make an coarser estimate on the database size; we can get a good estimate if
immediate input is specified (fasta data or a file). However, if a taxon list or
genome list is specfied this is tricker. For a genome list we just make an estimate based
on average genome sizes. For a taxon list, we do a lookup to find the blast database files
and compute based on that. (This size data should get pushed into data in C<database.json>.
Next time 'round.)

=cut

sub preflight
{
    my($self, $app, $app_def, $raw_params, $params) = @_;

    my($inp_type, $inp_size_est) = $self->compute_input_preflight($app, $params);
    my($db_type, $db_size_est) = $self->compute_db_preflight($app, $params);

    my $blast = $self->determine_blast_program($params);

    print "PF inp $inp_type $inp_size_est\n";
    print "PF db  $db_type $db_size_est\n";

    my $time;
    if ($db_size_est < 10_000)
    {
	$time = 1 * 60;
    }
    elsif ($db_size_est < 100_000)
    {
	$time = 5 * 60;
    }
    else
    {
	$time = 10 * 60;
    }
    my $inp_factor = ceil($inp_size_est / 1_000);

    $time *= $inp_factor;
    $time *= 6 if $blast eq 'tblastn' || $blast eq 'tblastx';

    #
    # Require at least 2 hours. If this run triggers the download of a new
    # container we will run out of time.
    #
    my $min_time = 120 * 60;
    
    $time = $min_time if $time < $min_time;

    #
    # Cap the runtime at 5 days. There are some glitches in long run compute
    # estimates that make the time estimate wildly wrong.
    #

    my $max_time = 86400 * 5;
    $time = $max_time if $time > $max_time;
    
    my $cpu = 4;

    $cpu = 8 if $db_size_est > 1_000_000;

    my $mem = "48G";
    if ($params->{db_source} eq 'taxon_list')
    {
	$mem = "128G";
    }

    my $pf = {
	cpu => $cpu,
	memory => "32G",
	runtime => $time,
	storage => 0,
	is_control_task => 0,
	policy_data => { constraint => 'sim' }
    };
    return $pf;
}

sub process
{
    my($self, $app, $app_def, $raw_params, $params) = @_;
    print "Proc Homology ", Dumper($app_def, $raw_params, $params);

    $self->app($app);
    $self->app_def($app_def);
    $self->params($params);
    $self->token($app->token);
    $self->task_id($app->task_id);

    my $cwd = getcwd();
    my $work_dir = "$cwd/work";
    my $stage_dir = "$cwd/stage";
    
    -d $work_dir or mkdir $work_dir or die "Cannot mkdir $work_dir: $!";
    -d $stage_dir or mkdir $stage_dir or die "Cannot mkdir $stage_dir: $!";

    $self->work_dir($work_dir);
    $self->stage_dir($stage_dir);

    my $output_base = $self->params->{output_file};
    my $output_folder = $self->app->result_folder();

    $self->output_base($output_base);
    $self->output_folder($output_folder);

    my($input_file, $is_short) = $self->stage_input($params, $stage_dir);
    my($db_file, $blast_params) = $self->stage_database($params, $stage_dir);

    if (!$db_file)
    {
	die "Unable to find database \n";
    }

    my $blast_program = $self->determine_blast_program($params);
    $self->blast_program($blast_program);

    $self->run_blast($params, $input_file, $db_file, $blast_program, $blast_params, $work_dir, $is_short);

    # Write output

    my %typemap = (json => 'json',
		   txt => 'txt',
		   archive => 'unspecified',
		   tsv => 'tsv');

    for my $out (glob("$work_dir/*"))
    {
	my $file = basename($out);
	my($suffix) = $out =~ /\.([^.]+)$/;
	my $type = $typemap{$suffix} // 'unspecified';
	$self->app->workspace->save_file_to_file($out, {}, "$output_folder/$file", $type, 1, 1);
    }
	
}

=head2 stage_input

    my $input_file = $self->stage_input($params, $stage_dir)

Stage the input as defined by C<$params> into C<$stage_dir>.

=cut

sub stage_input
{
    my($self, $params, $stage_dir) = @_;

    my $src = $params->{input_source};
    my $method = "stage_input_$src";
    my($file, $is_short) = $self->$method($params, $stage_dir);
    return($file, $is_short);
}

sub stage_input_id_list
{
    my($self, $params, $stage_dir) = @_;

    my $ids = $self->app->validate_param_array('input_id_list');
    if (!defined($ids))
    {
	die "Input source defined as id_list, but no input_id_list parameter specified";
    }

    return $self->stage_input_ids($params, $stage_dir, $ids, 0);
}

sub stage_input_feature_group
{
    my($self, $params, $stage_dir) = @_;

    my $group = $params->{input_feature_group};
    
    if (!defined($group))
    {
	die "Input source defined as feature_group, but no input_feature_group parameter specified";
    }

    my $ids = $self->api->retrieve_feature_ids_from_feature_group($group);
    printf STDERR "Loaded %d ids\n", scalar @$ids;

    return $self->stage_input_ids($params, $stage_dir, $ids, 1);
}

sub stage_input_ids
{
    my($self, $params, $stage_dir, $ids, $use_feature_ids) = @_;

    printf STDERR "Staging %d ids use_feature_ids='$use_feature_ids'\n", scalar @$ids;

    my $seqs;
    if ($params->{input_type} eq 'aa')
    {
	$seqs = $self->api->retrieve_protein_feature_sequence($ids, $use_feature_ids);
    }
    else
    {
	$seqs = $self->api->retrieve_nucleotide_feature_sequence($ids, $use_feature_ids);
    }
    printf STDERR "Staged %d sequences\n", scalar keys %$seqs;
    my $file = "$stage_dir/input_$params->{input_type}.fa";
    my $is_short = 1;
    if (open(my $fh, ">", $file))
    {
	#
	# If we came in with feature IDs, we can't use the original IDs since
	# we mapped them to patric IDs when creating the sequence data
	#
	my $idlist = $use_feature_ids ? [keys %$seqs] : $ids;
	for my $id (@$idlist)
	{
	    $is_short = 0 if length($seqs->{$id}) > $self->{short_feature_threshold};
	    print_alignment_as_fasta($fh, [$id, undef, $seqs->{$id}]);
	}
	close($fh);
    }
    else
    {
	die "Cannot open $file for writing: $!";
    }
    return ($file, $is_short);
}

sub stage_input_fasta_data
{
    my($self, $params, $stage_dir) = @_;

    my $file = "$stage_dir/input_$params->{input_type}.fa";
    open(my $infh, "<", \$params->{input_fasta_data})
	or die "Cannot open input fasta data for reading: $!";

    open(my $fh, ">", $file)
	or die "Cannot open $file for writing: $!";

    my $stats = $self->read_and_validate_fasta($infh, $params->{input_type}, $fh, 1);

    my $is_short = $stats->min() <= $self->{short_feature_threshold};
    return ($file, $is_short);
}

#
# Stage input from a workspace file
#
sub stage_input_fasta_file
{
    my($self, $params, $stage_dir) = @_;

    my $tmp = File::Temp->new;

    my $file = "$stage_dir/input_$params->{input_type}.fa";
    my $ws = Bio::P3::Workspace::WorkspaceClientExt->new;

    my $is_short = 1;

    $ws->copy_files_to_handles(1, undef, [[$params->{input_fasta_file}, $tmp]]);
    $tmp->seek(0, SEEK_SET);
    
    open(my $fh, ">", $file) or die "Cannot open $file for writing: $!";

    my $stats = $self->read_and_validate_fasta($tmp, $params->{input_type}, $fh, 1);

    my $is_short = $stats->min() <= $self->{short_feature_threshold};
    close($fh);

    return($file, $is_short);
}

=head2 stage_database

    ($db_file, $blast_params) $self->stage_database($params, $stage_dir)

Stage the database as defined by C<$params> into C<$stage_dir>.

We return the database file and an initial set of BLAST parameters. The BLAST parameters
may be required if the database specification in the input parameters includes
a subset of taxonomy IDs to be searched.

=cut

sub stage_database
{
    my($self, $params, $stage_dir) = @_;

    my $src = $params->{db_source};
    my $method = "stage_database_$src";

    my $blast_params = [];

    my $db_file = $self->$method($params, $stage_dir, $blast_params);
    return ($db_file, $blast_params);
}

=head2 stage_database_fasta_data

Stage the fasta data from C<$params->{db_fasta_data}> to a local file.

=cut

sub stage_database_fasta_data
{
    my($self, $params, $stage_dir, $blast_params) = @_;

    my $file = "$stage_dir/db_$params->{db_type}.fa";

    if (!exists $params->{db_fasta_data})
    {
	die "Database of type fasta_data requested but parameter db_fasta_data is missing\n";
    }

    open(my $fh, "<", \ $params->{db_fasta_data}) or die "Cannot open filehandle on data: $!";
    open(my $out, ">", $file) or die "Cannot open output file $file: $!";

    $self->read_and_validate_fasta($fh, $params->{db_type}, $out, 0);

    close($fh);
    close($out);

    my $ok = run(["makeblastdb",
		  "-dbtype", $makeblastdb_dbtype{$params->{db_type}},
		  "-in", $file]);
    $ok or die "makeblastdb failed: $!";

    return $file;
}

sub stage_database_fasta_file
{
    my($self, $params, $stage_dir, $blast_params) = @_;

    my $file = "$stage_dir/db_$params->{db_type}.fa";
    my $ws = Bio::P3::Workspace::WorkspaceClientExt->new;

    if (open(my $fh, ">", $file))
    {

	$ws->copy_files_to_handles(1, undef, [[$params->{db_fasta_file}, $fh]]);

	close($fh);
    }
    else
    {
	die "Cannot open $file for writing: $!";
    }

    my $ok = run(["makeblastdb",
		  "-dbtype", $makeblastdb_dbtype{$params->{db_type}},
		  "-in", $file]);
    $ok or die "makeblastdb failed: $!";

    return $file;
}

sub stage_database_ids
{
    my($self, $params, $stage_dir, $ids, $blast_params) = @_;

    my $seqs;
    if ($makeblastdb_dbtype{$params->{db_type}} eq 'prot')
    {
	$seqs = $self->api->retrieve_protein_feature_sequence($ids);
    }
    else
    {
	$seqs = $self->api->retrieve_nucleotide_feature_sequence($ids);
    }

    my $file = "$stage_dir/db_$params->{db_type}.fa";
    if (open(my $fh, ">", $file))
    {
	for my $id (@$ids)
	{
	    print_alignment_as_fasta($fh, [$id, undef, $seqs->{$id}]);
	}
	close($fh);
    }
    else
    {
	die "Cannot open $file for writing: $!";
    }
    my $ok = run(["makeblastdb",
		  "-dbtype", $makeblastdb_dbtype{$params->{db_type}},
		  "-in", $file]);
    $ok or die "makeblastdb failed: $!";
    return $file;
}

sub stage_database_id_list
{
    my($self, $params, $stage_dir, $blast_params) = @_;

    my $ids = $self->app->validate_param_array('db_id_list');
    if (!defined($ids))
    {
	die "Database source defined as id_list, but no db_id_list parameter specified";
    }

    my $file = $self->stage_database_ids($params, $stage_dir, $ids, $blast_params);

    return $file;
}

sub stage_database_feature_group
{
    my($self, $params, $stage_dir, $blast_params) = @_;

    my $group = $params->{db_feature_group};
    
    my $ids = $self->api->retrieve_patricids_from_feature_group($group);

    my $file = $self->stage_database_ids($params, $stage_dir, $ids, $blast_params);

    return $file;
}

sub stage_database_genome_list
{
    my($self, $params, $stage_dir, $blast_params) = @_;

    my $genomes = $params->{db_genome_list};

    $self->stage_database_genomes($params, $genomes, $stage_dir, $blast_params);
}

sub stage_database_genome_group
{
    my($self, $params, $stage_dir, $blast_params) = @_;

    my $group = $params->{db_genome_group};
    my $genomes = $self->api->retrieve_patric_ids_from_genome_group($group);
    
    $self->stage_database_genomes($params, $genomes, $stage_dir, $blast_params);
}

sub stage_database_genomes
{
    my($self, $params, $genomes, $stage_dir, $blast_params) = @_;

    my $dblist = "$stage_dir/genome-list";
    open(D, ">", $dblist) or die "Cannot write $dblist: $!";
    for my $g (@$genomes)
    {
	if ($g =~ /^(\d+\.\d+)$/)
	{
	    print D "$1\n";
	    print STDERR "Add genome: $1\n";
	}
	else
	{
	    die "Invalid genome id $g in genome list";
	}
    }
    close(D);
    if (! -s $dblist)
    {
	die "No genomes written from genome list";
    }

    my $db = "$stage_dir/genomes";
    my @cmd = ("p3-make-genome-list-blast-database",
	       "--ids-from", $dblist,
	       "--db-type", $params->{db_type},
	       $db);
    print STDERR "Run: @cmd\n";
    
    my $rc = system(@cmd);
    $rc == 0 or die "Error creating blast database";

    return $db;
}

=head2 stage_database_taxon_list

Handle a request for a list of taxa.

We search the database list for precomputed databases that contain the desired taxa.

The database list is structured as a list of object with two primary properties:

=over 4

=item name

The name of the database

=item db_list

A listing of the BLAST databases that comprise the database, with their properties.

=back

The C<db_list> is the list of available BLAST databases.

Each item has the following properties:

=over 4

=item path

The location in the filesystem of the BLAST database; relative a toplevel path
defined by the execution environment. 

=item type

The type of the database, as listed above under L<Database types>

=item ftype

Either "features" or "contigs"; this may be redundant.

=item genome_counts

An object mapping the genome id contained in the data to the count
of sequences with that genome  id.

=item tax_counts

An object mapping the taxon id contained in the data to the count
of sequences with that taxon id. 

=back

We find the database or databases to use by scanning the tax_counts fields in the available
database of the correct type for the presence of the desired taxa. 

=cut
    
sub stage_database_taxon_list
{
    my($self, $params, $stage_dir, $blast_params) = @_;

    my $dbtype = $params->{db_type};

    my $db_list = $params->{db_taxon_list};
    if (!$db_list || ref($db_list) ne 'ARRAY')
    {
	die "stage_database_taxon_list: invalid or missing db_taxon_list";
    }

    my($taxa, $file_list) = $self->blastdbs->find_databases_for_taxa($dbtype, $db_list);

    if (!$taxa)
    {
	die "No taxa were found to search";
    }

    if (@$taxa > 10)
    {
	my $tlist = $self->stage_dir . "/taxids";

	open(T, ">", $tlist) or die "Cannot write $tlist: $!";
	print T "$_\n" foreach @$taxa;
	close(T);
	
	push(@$blast_params, "-taxidlist", $tlist);
    }
    else
    {
	push(@$blast_params, "-taxids", join(",", @$taxa));
    }

    #
    # If there are multiples, construct an alias database.
    #
    
    if (@$file_list > 1)
    {
	my $alias = $self->stage_dir . "/aliasdb";
	my $files = $self->stage_dir . "/dbfiles";

	open(my $fh, ">", $files) or die "Cannot write $files: $!";
	print $fh "$_\n" foreach @$file_list;
	close($fh);
	
	my $ok = run(["blastdb_aliastool",
		      "-dblist_file", $files,
		      "-dbtype", $makeblastdb_dbtype{$dbtype},
		      "-title", "combined taxon database",
		      "-out", $alias]);
	if (!$ok)
	{
	    die "Cannot create alias database for desired taxa @{$params->{db_taxon_list}}\n";
	}

	return $alias;
    }
    else
    {
	return $file_list->[0];
    }
}

sub stage_database_precomputed_database
{
    my($self, $params, $stage_dir, $blast_params) = @_;

    my $db = $params->{db_precomputed_database};
    defined($db) or die "Precomputed database was selected but db_precomputed_database key was missing from parameters";

    return $self->blastdbs->find_precomputed_database($db, $params->{db_type});
}

sub read_and_validate_fasta
{
    my($self, $in_fh, $type, $out_fh, $collect_stats) = @_;

    my $is_prot = $type eq 'faa' || $type eq 'aa';

    my $stats;
    if ($collect_stats)
    {
	$stats = Statistics::Descriptive::Sparse->new();
    }

    parse_fasta($in_fh, undef, sub {
	my($id, $seq) = @_;
	$stats->add_data(length($seq)) if $stats;
	print_alignment_as_fasta($out_fh, [$id, undef, $seq]);
	return 1;
    }, $is_prot);
    return $stats;
}

=head2 determine_blast_program
    
Determine the appropriate BLAST program based on the input and db type and
a setting for C<blast_program> from the parameters if present.
    
=cut

sub determine_blast_program
{
    my($self, $params) = @_;
    
    my $valid_programs = $search_type_map{$params->{input_type}}->{$params->{db_type}};
    if (!$valid_programs)
    {
	die "No valid program found for input type '$params->{input_type}' and db type '$params->{db_type}'\n";
    }
    my $prog = $params->{blast_program};
    if ($prog)
    {
	if (none { $_ eq $prog } @$valid_programs)
	{
	    die "$prog is not a valid blast program for input type '$params->{input_type}' and db type '$params->{db_type}'\n";
	}
    }
    else
    {
	$prog = $valid_programs->[0];
    }
    return $prog;
}

sub run_blast
{
    my($self, $params, $input_file, $db_file, $blast_program, $blast_params, $work_dir, $is_short) = @_;

    #
    # Set output format to blast archive, then translate to desired formats.
    #

    my @opts = @$blast_params;

    my $cpus = $ENV{P3_ALLOCATED_CPU} // 1;

    push(@opts, "-num_threads", $cpus);

    my $out_file = "$work_dir/blast_out.archive";
    my $out_json_file = "$work_dir/blast_out.raw.json";
    my $out_tbl_file = "$work_dir/blast_out.txt";
    my $out_tbl_hdrs = "$work_dir/blast_headers.txt";
    my $out_proc_json_file = "$work_dir/blast_out.json";
    my $out_md_file = "$work_dir/blast_out.metadata.json";

    push(@opts,
	 "-query", $input_file,
	 "-db", $db_file,
	 "-outfmt", 11,
	 "-out", $out_file,
	);

    if ($blast_program eq 'blastn' && $is_short)
    {
	push(@opts, "-task", "blastn-short");
    }
    elsif ($blast_program eq 'blastp' && $is_short)
    {
	push(@opts, "-task", "blastp-short");
    }


    if (exists $params->{blast_evalue_cutoff})
    {
	push(@opts, "-evalue", $params->{blast_evalue_cutoff});
    }

    if (exists $params->{blast_max_hits})
    {
	push(@opts, "-max_target_seqs", $params->{blast_max_hits});
    }

    if (exists $params->{blast_min_coverage})
    {
	push(@opts, "-qcov_hsp_perc", $params->{blast_min_coverage});
    }

    print  "$blast_program @opts\n";
    my $ok = run([$blast_program, @opts]);
    $ok or die "Error running $blast_program @opts: $!";

    #
    # We wrote blast archive output. Use blast_formatter to
    # create json output.
    #

    my $ok = run(["blast_formatter",
	       "-archive", $out_file,
	       "-outfmt", 15,
	       "-out", $out_json_file]);
    $ok or die "Error running blast formatter: $!";

    #
    # Create tabular output; we need to filter to remove gnl| from the ids.
    #
    my @cols = qw(qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen);

    write_file($out_tbl_hdrs, join("\t", @cols) . "\n");

    my $tbl_fh;
    my $pipe = IO::Handle->new;
    
    open($tbl_fh, ">", $out_tbl_file) or die "Cannot write $out_tbl_file: $!";
    my $h = start(["blast_formatter",
		   "-archive", $out_file,
		   "-outfmt", join(" ", 6, @cols)],
		  ">pipe", $pipe);
    $h or die "Error running blast formatter: $!";

    while (<$pipe>)
    {
	s/^([^\t]+\t)gnl\|/$1/;
	print $tbl_fh $_;
    }

    finish($h) or die "returned: $?";
    
    close($pipe);
    close($tbl_fh);

    my $txt = read_file($out_json_file);
    my($doc, $metadata) = $self->massage_blast_json($txt);
    write_file($out_proc_json_file, $self->json->encode($doc));
    write_file($out_md_file, $self->json->encode($metadata));
}

sub massage_blast_json
{
    my($self, $json) = @_;
    my $doc = eval { $self->json->decode($json) };
    if ($@)
    {
	die "json parse failed: $@\n$json\n";
    }
    $doc = $doc->{BlastOutput2};
    $doc or die "JSON output didn't have expected key BlastOutput2";

    my $metadata = {};
    for my $report (@$doc)
    {
	my $search = $report->{report}->{results}->{search};
	$search->{query_id} =~ s/^gnl\|//;
	if ($search->{query_id} =~ /^Query_\d+/)
	{
	    my($xid) = $search->{query_title} =~ /^(\S+)/;
	    $search->{query_id} = $xid if $xid;
	}
	for my $res (@{$search->{hits}})
	{
	    for my $desc (@{$res->{description}})
	    {
		my $md = $self->decode_title($desc);
		
		$metadata->{$desc->{id}} = $md if $md;
	    }
	}
    }
    return($doc, $metadata);
}

sub decode_title
{
    my($self, $desc) = @_;
	
    my $md;
    
    if ($desc->{id} =~ /^gnl\|BL_ORD/)
    {
	if ($desc->{title} =~ /^(\S+)\s{3}(.*)\s{3}(.*)$/)
	{
	    my $id = $1;
	    my $fun = $2;
	    my $ginfo = $3;
	    $id =~ s/^((fig\|\d+\.\d+.[^.]+\.\d+)|(accn\|[^|]+))//;
	    my $fid = $1;
	    print "id=$id\n";
	    $id =~ s/^\|//;
	    $id =~ s/\|$//;
	    my @rest = split(/\|/, $id);
	    my $locus;
	    my $alt;
	    if (@rest == 2)
	    {
		($locus, $alt) = @rest;
		$md->{locus_tag} = $locus;
		$md->{alt_locus_tag} = $alt;
	    }
	    elsif (@rest == 1)
	    {
		$alt = $rest[0];
		$md->{alt_locus_tag} = $alt;
	    }
	    if ($ginfo =~ /^\[(.*) \| (\d+\.\d+)/)
	    {
		$md->{genome_name} = $1;
		$md->{genome_id} = $2;
	    }
	    $desc->{id} = $fid;
	    $md->{function} = $fun;
	}
	elsif ($desc->{title} =~ /^(\S+)\s+(.*)\s{3}\[(.*?)(\s*\|\s*(\S+))?\]\s*$/)
	{
	    $desc->{id} = $1;
	    $md->{function} = $2;
	    $md->{genome_name} = $3;
	    $md->{genome_id} = $5 if $5;
	}
	elsif ($desc->{title} =~ /^(\S+)\s+\[(.*?)(\s*\|\s*(\S+))?\]\s*$/)
	{
	    $desc->{id} = $1;
	    $md->{genome_name} = $2;
	    $md->{genome_id} = $4 if $4;
	}
	elsif ($desc->{title} =~ /^(\S+)\s+(.*)\s*$/)
	{
	    $desc->{id} = $1;
	    $md->{function} = $2;
	    $desc->{title} = $2;
	}
	else
	{
	    $desc->{id} = $desc->{title};
	}
    }
    else
    {
	$desc->{id} =~ s/^gnl\|//;
	if ($desc->{title} =~ /^\s*(.*)\s{3}\[(.*?)(\s*\|\s*(\S+))?\]\s*$/)
	{
	    $md->{function} = $1;
	    $md->{genome_name} = $2;
	    $md->{genome_id} = $4 if $4;
	}
	elsif ($desc->{title} =~ /^\s*\[(.*?)(\s*\|\s*(\S+))?\]\s*$/)
	{
	    $md->{genome_name} = $1;
	    $md->{genome_id} = $3 if $3;
	}
    }
    if ($desc->{id} =~ /^(kb\|g\.\d+)/)
    {
	$md->{genome_id} = $1;
    }
    return $md;
}

sub find_db_in_path
{
    my($self, $path) = @_;
    
    for my $p (@{$self->{database_path}})
    {
	my $file = "$p/$path";
	for my $suf (qw(pal nal psq nsq))
	{
	    my $chk = "$file.$suf";
	    if (-f $chk)
	    {
		return ($file);
	    }
	}
    }
    
}

sub compute_input_preflight
{
    my($self, $app, $params) = @_;
    # ret   my($inp_type, $inp_size_est)

    my $inp_type = $params->{input_type};

    my $src = $params->{input_source};
    my $sz;
    if ($src eq 'id_list')
    {
	$sz = @{$params->{input_id_list}} * 1000;
	$sz *= 3 if $inp_type eq 'dna';
    }
    elsif ($src eq 'fasta_data')
    {
	$sz = length($params->{input_fasta_data});
    }
    elsif ($src eq 'fasta_file')
    {
	my $ws = $app->workspace();
	my $stat = $ws->stat($params->{input_fasta_file});
	if ($stat)
	{
	    $sz = $stat->size;
	}
	else
	{
	    die "Input file $params->{input_fasta_file} not found\n";
	}
    }
    elsif ($src eq 'feature_group')
    {
	my $ids = $self->api->retrieve_feature_ids_from_feature_group($params->{input_feature_group});
	$sz = @$ids * 1000;
	$sz *= 3 if $inp_type eq 'dna';
    }
    else
    {
	die "Invalid input type $src\n";
    }

    return($inp_type, $sz);

}


sub compute_db_preflight
{
    my($self, $app, $params) = @_;

    my $db_type = $params->{db_type};
    my $db_src  = $params->{db_source};

    my $sz;

    if ($db_src eq 'fasta_data')
    {
	$sz = length($params->{db_fasta_data});
    }
    elsif ($db_src eq 'fasta_file')
    {
	my $ws = $app->workspace();
	my $stat = $ws->stat($params->{db_fasta_file});
	if ($stat)
	{
	    $sz = $stat->size;
	}
	else
	{
	    die "Database file $params->{db_fasta_file} not found\n";
	}
    }
    elsif ($db_src eq 'genome_list')
    {
	$sz = 1000000;
	$sz *= 3 if $db_type ne 'faa';
    }
    elsif ($db_src eq 'genome_group')
    {
	$sz = 1000000;
	$sz *= 3 if $db_type ne 'faa';
    }
    elsif ($db_src eq 'feature_group')
    {
	$sz = 1000000;
	$sz *= 3 if $db_type ne 'faa';
    }
    elsif ($db_src eq 'taxon_list')
    {
	my($taxa, $file_list) = $self->blastdbs->find_databases_for_taxa($db_type, $params->{db_taxon_list});
	for my $f (@$file_list)
	{
	    for my $bf (glob("$f*[pn]sq"))
	    {
		$sz += -s $bf;
	    }
	}
    }
    elsif ($db_src eq 'precomputed_database')
    {
	my $f = $self->blastdbs->find_precomputed_database($params->{db_precomputed_database}, $db_type);
	for my $bf (glob("$f*[pn]sq"))
	{
	    my $fsize = -s $bf;
	    $sz += $fsize;
	}
    }
    else
    {
	die "Unknown database source $db_src\n";
    }

    return($db_type, $sz);
}
1;
