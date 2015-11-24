package PangenomeOrthomcl::PangenomeOrthomclClient;

use JSON::RPC::Client;
use POSIX;
use strict;
use Data::Dumper;
use URI;
use Bio::KBase::Exceptions;
my $get_time = sub { time, 0 };
eval {
    require Time::HiRes;
    $get_time = sub { Time::HiRes::gettimeofday() };
};

use Bio::KBase::AuthToken;

# Client version should match Impl version
# This is a Semantic Version number,
# http://semver.org
our $VERSION = "0.1.0";

=head1 NAME

PangenomeOrthomcl::PangenomeOrthomclClient

=head1 DESCRIPTION


A KBase module: PangenomeOrthomcl


=cut

sub new
{
    my($class, $url, @args) = @_;
    

    my $self = {
	client => PangenomeOrthomcl::PangenomeOrthomclClient::RpcClient->new,
	url => $url,
	headers => [],
    };

    chomp($self->{hostname} = `hostname`);
    $self->{hostname} ||= 'unknown-host';

    #
    # Set up for propagating KBRPC_TAG and KBRPC_METADATA environment variables through
    # to invoked services. If these values are not set, we create a new tag
    # and a metadata field with basic information about the invoking script.
    #
    if ($ENV{KBRPC_TAG})
    {
	$self->{kbrpc_tag} = $ENV{KBRPC_TAG};
    }
    else
    {
	my ($t, $us) = &$get_time();
	$us = sprintf("%06d", $us);
	my $ts = strftime("%Y-%m-%dT%H:%M:%S.${us}Z", gmtime $t);
	$self->{kbrpc_tag} = "C:$0:$self->{hostname}:$$:$ts";
    }
    push(@{$self->{headers}}, 'Kbrpc-Tag', $self->{kbrpc_tag});

    if ($ENV{KBRPC_METADATA})
    {
	$self->{kbrpc_metadata} = $ENV{KBRPC_METADATA};
	push(@{$self->{headers}}, 'Kbrpc-Metadata', $self->{kbrpc_metadata});
    }

    if ($ENV{KBRPC_ERROR_DEST})
    {
	$self->{kbrpc_error_dest} = $ENV{KBRPC_ERROR_DEST};
	push(@{$self->{headers}}, 'Kbrpc-Errordest', $self->{kbrpc_error_dest});
    }

    #
    # This module requires authentication.
    #
    # We create an auth token, passing through the arguments that we were (hopefully) given.

    {
	my $token = Bio::KBase::AuthToken->new(@args);
	
	if (!$token->error_message)
	{
	    $self->{token} = $token->token;
	    $self->{client}->{token} = $token->token;
	}
        else
        {
	    #
	    # All methods in this module require authentication. In this case, if we
	    # don't have a token, we can't continue.
	    #
	    die "Authentication failed: " . $token->error_message;
	}
    }

    my $ua = $self->{client}->ua;	 
    my $timeout = $ENV{CDMI_TIMEOUT} || (30 * 60);	 
    $ua->timeout($timeout);
    bless $self, $class;
    #    $self->_validate_version();
    return $self;
}




=head2 build_pangenome_with_orthomcl

  $return = $obj->build_pangenome_with_orthomcl($params)

=over 4

=item Parameter and return types

=begin html

<pre>
$params is a PangenomeOrthomcl.BuildPangenomeWithOrthmclParams
$return is a PangenomeOrthomcl.BuildPangenomeWithOrthmclResult
BuildPangenomeWithOrthmclParams is a reference to a hash where the following keys are defined:
	input_genomeset_ref has a value which is a PangenomeOrthomcl.ws_genomeset_id
	input_genome_refs has a value which is a reference to a list where each element is a PangenomeOrthomcl.ws_genome_id
	output_workspace has a value which is a PangenomeOrthomcl.workspace
	output_pangenome_id has a value which is a PangenomeOrthomcl.ws_pangenome_id
	num_descriptions has a value which is an int
	num_alignments has a value which is an int
	evalue has a value which is a string
	word_size has a value which is an int
	gapopen has a value which is an int
	gapextend has a value which is an int
	matrix has a value which is a string
	threshold has a value which is an int
	comp_based_stats has a value which is a string
	seg has a value which is a string
	lcase_masking has a value which is a PangenomeOrthomcl.boolean
	xdrop_gap_final has a value which is a float
	window_size has a value which is an int
	use_sw_tback has a value which is a PangenomeOrthomcl.boolean
	mcl_p has a value which is an int
	mcl_s has a value which is an int
	mcl_r has a value which is an int
	mcl_pct has a value which is an int
	mcl_warn_p has a value which is an int
	mcl_warn_factor has a value which is an int
	mcl_init_l has a value which is an int
	mcl_main_l has a value which is an int
	mcl_init_i has a value which is a float
	mcl_main_i has a value which is a float
ws_genomeset_id is a string
ws_genome_id is a string
workspace is a string
ws_pangenome_id is a string
boolean is an int
BuildPangenomeWithOrthmclResult is a reference to a hash where the following keys are defined:
	output_log has a value which is a string
	pangenome_ref has a value which is a PangenomeOrthomcl.ws_pangenome_id

</pre>

=end html

=begin text

$params is a PangenomeOrthomcl.BuildPangenomeWithOrthmclParams
$return is a PangenomeOrthomcl.BuildPangenomeWithOrthmclResult
BuildPangenomeWithOrthmclParams is a reference to a hash where the following keys are defined:
	input_genomeset_ref has a value which is a PangenomeOrthomcl.ws_genomeset_id
	input_genome_refs has a value which is a reference to a list where each element is a PangenomeOrthomcl.ws_genome_id
	output_workspace has a value which is a PangenomeOrthomcl.workspace
	output_pangenome_id has a value which is a PangenomeOrthomcl.ws_pangenome_id
	num_descriptions has a value which is an int
	num_alignments has a value which is an int
	evalue has a value which is a string
	word_size has a value which is an int
	gapopen has a value which is an int
	gapextend has a value which is an int
	matrix has a value which is a string
	threshold has a value which is an int
	comp_based_stats has a value which is a string
	seg has a value which is a string
	lcase_masking has a value which is a PangenomeOrthomcl.boolean
	xdrop_gap_final has a value which is a float
	window_size has a value which is an int
	use_sw_tback has a value which is a PangenomeOrthomcl.boolean
	mcl_p has a value which is an int
	mcl_s has a value which is an int
	mcl_r has a value which is an int
	mcl_pct has a value which is an int
	mcl_warn_p has a value which is an int
	mcl_warn_factor has a value which is an int
	mcl_init_l has a value which is an int
	mcl_main_l has a value which is an int
	mcl_init_i has a value which is a float
	mcl_main_i has a value which is a float
ws_genomeset_id is a string
ws_genome_id is a string
workspace is a string
ws_pangenome_id is a string
boolean is an int
BuildPangenomeWithOrthmclResult is a reference to a hash where the following keys are defined:
	output_log has a value which is a string
	pangenome_ref has a value which is a PangenomeOrthomcl.ws_pangenome_id


=end text

=item Description



=back

=cut

 sub build_pangenome_with_orthomcl
{
    my($self, @args) = @_;

# Authentication: required

    if ((my $n = @args) != 1)
    {
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error =>
							       "Invalid argument count for function build_pangenome_with_orthomcl (received $n, expecting 1)");
    }
    {
	my($params) = @args;

	my @_bad_arguments;
        (ref($params) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument 1 \"params\" (value was \"$params\")");
        if (@_bad_arguments) {
	    my $msg = "Invalid arguments passed to build_pangenome_with_orthomcl:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	    Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
								   method_name => 'build_pangenome_with_orthomcl');
	}
    }

    my $result = $self->{client}->call($self->{url}, $self->{headers}, {
	method => "PangenomeOrthomcl.build_pangenome_with_orthomcl",
	params => \@args,
    });
    if ($result) {
	if ($result->is_error) {
	    Bio::KBase::Exceptions::JSONRPC->throw(error => $result->error_message,
					       code => $result->content->{error}->{code},
					       method_name => 'build_pangenome_with_orthomcl',
					       data => $result->content->{error}->{error} # JSON::RPC::ReturnObject only supports JSONRPC 1.1 or 1.O
					      );
	} else {
	    return wantarray ? @{$result->result} : $result->result->[0];
	}
    } else {
        Bio::KBase::Exceptions::HTTP->throw(error => "Error invoking method build_pangenome_with_orthomcl",
					    status_line => $self->{client}->status_line,
					    method_name => 'build_pangenome_with_orthomcl',
				       );
    }
}
 
  

sub version {
    my ($self) = @_;
    my $result = $self->{client}->call($self->{url}, $self->{headers}, {
        method => "PangenomeOrthomcl.version",
        params => [],
    });
    if ($result) {
        if ($result->is_error) {
            Bio::KBase::Exceptions::JSONRPC->throw(
                error => $result->error_message,
                code => $result->content->{code},
                method_name => 'build_pangenome_with_orthomcl',
            );
        } else {
            return wantarray ? @{$result->result} : $result->result->[0];
        }
    } else {
        Bio::KBase::Exceptions::HTTP->throw(
            error => "Error invoking method build_pangenome_with_orthomcl",
            status_line => $self->{client}->status_line,
            method_name => 'build_pangenome_with_orthomcl',
        );
    }
}

sub _validate_version {
    my ($self) = @_;
    my $svr_version = $self->version();
    my $client_version = $VERSION;
    my ($cMajor, $cMinor) = split(/\./, $client_version);
    my ($sMajor, $sMinor) = split(/\./, $svr_version);
    if ($sMajor != $cMajor) {
        Bio::KBase::Exceptions::ClientServerIncompatible->throw(
            error => "Major version numbers differ.",
            server_version => $svr_version,
            client_version => $client_version
        );
    }
    if ($sMinor < $cMinor) {
        Bio::KBase::Exceptions::ClientServerIncompatible->throw(
            error => "Client minor version greater than Server minor version.",
            server_version => $svr_version,
            client_version => $client_version
        );
    }
    if ($sMinor > $cMinor) {
        warn "New client version available for PangenomeOrthomcl::PangenomeOrthomclClient\n";
    }
    if ($sMajor == 0) {
        warn "PangenomeOrthomcl::PangenomeOrthomclClient version is $svr_version. API subject to change.\n";
    }
}

=head1 TYPES



=head2 boolean

=over 4



=item Description

Indicates true or false values, false = 0, true = 1
@range [0,1]


=item Definition

=begin html

<pre>
an int
</pre>

=end html

=begin text

an int

=end text

=back



=head2 workspace

=over 4



=item Description

Name of workspace.


=item Definition

=begin html

<pre>
a string
</pre>

=end html

=begin text

a string

=end text

=back



=head2 ws_genome_id

=over 4



=item Description

The workspace ID for a GenomeSet data object.
@id ws KBaseGenome.Genome


=item Definition

=begin html

<pre>
a string
</pre>

=end html

=begin text

a string

=end text

=back



=head2 ws_genomeset_id

=over 4



=item Description

The workspace ID for a GenomeSet data object.
@id ws KBaseSearch.GenomeSet


=item Definition

=begin html

<pre>
a string
</pre>

=end html

=begin text

a string

=end text

=back



=head2 ws_pangenome_id

=over 4



=item Description

The workspace ID for a Pangenome data object.
@id ws KBaseGenomes.Pangenome


=item Definition

=begin html

<pre>
a string
</pre>

=end html

=begin text

a string

=end text

=back



=head2 BuildPangenomeWithOrthmclParams

=over 4



=item Description

Input parameters of build_pangenome_with_orthomcl method.
input_genomeset_ref - optional input reference to genome set 
    object (alternative way is input_genome_refs);
input_genome_refs - optional input list of references to
    genome objects (alternative way is input_genomeset_ref);
output_workspace - workspace for saving resulting pangenome;
output_pangenome_id - name of resulting pangenome object;
num_descriptions - [blastp, -v] Store one-line descriptions for 
    this number of database sequences. Default value is 100000.
num_alignments - [blastp, -b] Store alignments for this number of 
    database sequences. Default value is 100000.
evalue - [blastp, -e] Expect value (E) for saving hits. Default
    value is 1e-5.
word_size - [blastp, -W] Word size of initial match. Valid word 
    sizes are 2-7. Default value is 3.
gapopen - [blastp, -G] Cost to open a gap. Default value is 11.
gapextend - [blastp, -E] Cost to extend a gap. Default value is 1.
matrix - [blastp, -M] Scoring matrix name. Default value is BLOSUM62.
threshold - [blastp, -f] Minimum score to add a word to the BLAST 
    lookup table. Default value is 11.
comp_based_stats - [blastp, -C] Use composition-based statistics 
    (0: no composition-based statistics; 1: Composition-based 
    statistics as in NAR 29:2994-3005, 2001; 2: Composition-based 
    score adjustments as in Bioinformatics 21:902-911, 2005, 
    conditioned on sequence properties; 3 - Composition-based 
    score adjustment as in Bioinformatics 21:902-911, 2005, 
    unconditionally). Default value is 2.
seg - [blastp, -F] Filter query sequence with SEG (yes/no). Default
    value is yes.
lcase_masking - [blastp, -U] Use lower case filtering in query and 
    subject sequence(s). Default value is false(0).
xdrop_gap_final - [blastp, -Z] Heuristic value (in bits) for final 
    gapped alignment. Default value is 25.
window_size - [blastp, -A] Multiple hits window size, use 0 to 
    specify 1-hit algorithm. Default value is 40.
use_sw_tback - [blastp, -s] Compute locally optimal Smith-Waterman 
    alignments. Default value is false(0).
mcl_p - [mcl, -P] Prune number. Default value is 10000.
mcl_s - [mcl, -S] Selection number. Default value is 1100.
mcl_r - [mcl, -R] Recovery number. Default value is 1400.
mcl_pct - [mcl, -pct] Recovery percentage. Default value is 90.
mcl_warn_p - [mcl, -warn-p] Warn if pruning reduces mass to this 
    weight. Default value is 10.
mcl_warn_factor - [mcl, -warn-factor] Warn if pruning reduces entry 
    count by this value. Default value is 1000.
mcl_init_l - [mcl, -l] Initial loop length. Default value is 0.
mcl_main_l - [mcl, -L] Main loop length. Default value is 10000.
mcl_init_i - [mcl, -i] Initial inflation. Default value is 2.0.
mcl_main_i - [mcl, -I] Main inflation. Default value is 1.5.


=item Definition

=begin html

<pre>
a reference to a hash where the following keys are defined:
input_genomeset_ref has a value which is a PangenomeOrthomcl.ws_genomeset_id
input_genome_refs has a value which is a reference to a list where each element is a PangenomeOrthomcl.ws_genome_id
output_workspace has a value which is a PangenomeOrthomcl.workspace
output_pangenome_id has a value which is a PangenomeOrthomcl.ws_pangenome_id
num_descriptions has a value which is an int
num_alignments has a value which is an int
evalue has a value which is a string
word_size has a value which is an int
gapopen has a value which is an int
gapextend has a value which is an int
matrix has a value which is a string
threshold has a value which is an int
comp_based_stats has a value which is a string
seg has a value which is a string
lcase_masking has a value which is a PangenomeOrthomcl.boolean
xdrop_gap_final has a value which is a float
window_size has a value which is an int
use_sw_tback has a value which is a PangenomeOrthomcl.boolean
mcl_p has a value which is an int
mcl_s has a value which is an int
mcl_r has a value which is an int
mcl_pct has a value which is an int
mcl_warn_p has a value which is an int
mcl_warn_factor has a value which is an int
mcl_init_l has a value which is an int
mcl_main_l has a value which is an int
mcl_init_i has a value which is a float
mcl_main_i has a value which is a float

</pre>

=end html

=begin text

a reference to a hash where the following keys are defined:
input_genomeset_ref has a value which is a PangenomeOrthomcl.ws_genomeset_id
input_genome_refs has a value which is a reference to a list where each element is a PangenomeOrthomcl.ws_genome_id
output_workspace has a value which is a PangenomeOrthomcl.workspace
output_pangenome_id has a value which is a PangenomeOrthomcl.ws_pangenome_id
num_descriptions has a value which is an int
num_alignments has a value which is an int
evalue has a value which is a string
word_size has a value which is an int
gapopen has a value which is an int
gapextend has a value which is an int
matrix has a value which is a string
threshold has a value which is an int
comp_based_stats has a value which is a string
seg has a value which is a string
lcase_masking has a value which is a PangenomeOrthomcl.boolean
xdrop_gap_final has a value which is a float
window_size has a value which is an int
use_sw_tback has a value which is a PangenomeOrthomcl.boolean
mcl_p has a value which is an int
mcl_s has a value which is an int
mcl_r has a value which is an int
mcl_pct has a value which is an int
mcl_warn_p has a value which is an int
mcl_warn_factor has a value which is an int
mcl_init_l has a value which is an int
mcl_main_l has a value which is an int
mcl_init_i has a value which is a float
mcl_main_i has a value which is a float


=end text

=back



=head2 BuildPangenomeWithOrthmclResult

=over 4



=item Description

Output results of build_pangenome_with_orthomcl method.
One of 'pangenome_ref' and 'error' fields should be defined.


=item Definition

=begin html

<pre>
a reference to a hash where the following keys are defined:
output_log has a value which is a string
pangenome_ref has a value which is a PangenomeOrthomcl.ws_pangenome_id

</pre>

=end html

=begin text

a reference to a hash where the following keys are defined:
output_log has a value which is a string
pangenome_ref has a value which is a PangenomeOrthomcl.ws_pangenome_id


=end text

=back



=cut

package PangenomeOrthomcl::PangenomeOrthomclClient::RpcClient;
use base 'JSON::RPC::Client';
use POSIX;
use strict;

#
# Override JSON::RPC::Client::call because it doesn't handle error returns properly.
#

sub call {
    my ($self, $uri, $headers, $obj) = @_;
    my $result;


    {
	if ($uri =~ /\?/) {
	    $result = $self->_get($uri);
	}
	else {
	    Carp::croak "not hashref." unless (ref $obj eq 'HASH');
	    $result = $self->_post($uri, $headers, $obj);
	}

    }

    my $service = $obj->{method} =~ /^system\./ if ( $obj );

    $self->status_line($result->status_line);

    if ($result->is_success) {

        return unless($result->content); # notification?

        if ($service) {
            return JSON::RPC::ServiceObject->new($result, $self->json);
        }

        return JSON::RPC::ReturnObject->new($result, $self->json);
    }
    elsif ($result->content_type eq 'application/json')
    {
        return JSON::RPC::ReturnObject->new($result, $self->json);
    }
    else {
        return;
    }
}


sub _post {
    my ($self, $uri, $headers, $obj) = @_;
    my $json = $self->json;

    $obj->{version} ||= $self->{version} || '1.1';

    if ($obj->{version} eq '1.0') {
        delete $obj->{version};
        if (exists $obj->{id}) {
            $self->id($obj->{id}) if ($obj->{id}); # if undef, it is notification.
        }
        else {
            $obj->{id} = $self->id || ($self->id('JSON::RPC::Client'));
        }
    }
    else {
        # $obj->{id} = $self->id if (defined $self->id);
	# Assign a random number to the id if one hasn't been set
	$obj->{id} = (defined $self->id) ? $self->id : substr(rand(),2);
    }

    my $content = $json->encode($obj);

    $self->ua->post(
        $uri,
        Content_Type   => $self->{content_type},
        Content        => $content,
        Accept         => 'application/json',
	@$headers,
	($self->{token} ? (Authorization => $self->{token}) : ()),
    );
}



1;
