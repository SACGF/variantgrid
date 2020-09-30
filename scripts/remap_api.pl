# Uses modules:
#
# http://search.cpan.org/~jv/Getopt-Long-2.42/lib/Getopt/Long.pm
# http://search.cpan.org/~gaas/libwww-perl-6.05/lib/LWP/UserAgent.pm
# http://search.cpan.org/~gaas/HTTP-Message-6.06/lib/HTTP/Request/Common.pm
# http://search.cpan.org/~gaas/HTTP-Message-6.06/lib/HTTP/Headers.pm#
# http://search.cpan.org/~msergeant/XML-XPath-1.13/XPath.pm
# http://search.cpan.org/~msergeant/XML-XPath-1.13/XPath/XMLParser.pm
# http://search.cpan.org/~makamaka/JSON-2.90/lib/JSON.pm

use strict;
use warnings;

#for (@INC) { undef $_ }

require Getopt::Long;
require LWP::UserAgent;
require LWP::Protocol::https;
require HTTP::Request::Common;
require HTTP::Headers;
require XML::XPath;
require XML::XPath::XMLParser;
require JSON;


use Getopt::Long;
use LWP::UserAgent;
use LWP::Protocol::https;
use HTTP::Request::Common qw(POST), qw(GET);
use HTTP::Headers;
use XML::XPath;
use XML::XPath::XMLParser;
use JSON;


#use Data::Dumper;


my $protocol = "https";
my $netloc = "www.ncbi.nlm.nih.gov:443";

#my $URLbase = "http://".$netloc."/projects/genome/assembly/grc/TpfSub/";
my $remapdir = "/projects/genome/tools/remap/";
#my $remapdir = "/genome/tools/remap/";
my $URLbase = $protocol."://".$netloc."/".$remapdir;
my $RemapServiceCgi = "remapservice.cgi";
my $RemapInfoCgi = "remapinfo.cgi";

my $forward_url = "NOT_YET_SET";
my $jobid = "";
my $ctgtime = "";

my $annot_url = "";
my $report_url = "";
my $gbench_url = "";

#foreach (@INC) {
#    print `find $_ -type f -name "*.pm"`;
#}

# Install to     /am/ftp-pub-remap/remap_api.pl
# Download from  ftp://ftp.ncbi.nlm.nih.gov/pub/remap/remap_api.pl


# Usage: perl remap_api.pl  
#
#  Please limit your submissions to files of approximatly 250,000 rows, 
#  and at most 4 simultaneous submissions.
#
#  --mode <mode>
#    batches, asm-asm, asm-rsg, rsg-asm.
#  
#    batches will print a list of all assembly-assembly pairs that alignments 
#    exist for.
#
#  batches use no other argumentes, other modes
#  
#  required:
#  --annotation <filename> 
#    The file containing your existing annotation (GFF3, GVF, BED, etc)
#
#  --from <gencoll accession>
#  --dest <gencoll accession>
#    The gencoll accessions for the assemblies being mapped to and from
#    Use 'RefSeqGene' when mapping to RefSeqGene
#
#  options:
#  --allowdupes <'on' or 'off'>
#    on (default) does all possible mappings,
#    off excludes second pass mappings
#
#  --merge <'on' or 'off'>
#    on (default) will merge fragmented mappings back together
#    off will leave fragmented results untouched
#
#  --mincov <fraction 0.01 to 10.0, default 0.5 >
#  --maxexp <fraction 0.01 to 10.0, default 2.0 >
#    (r_cov = mapped_feature_length / original_feature_length 
#     if r_cov < min_cov then throw out feature
#     if r_cov > max_exp then throw out feature )
#    mapped_feature_length might not be identical to original_feature_length
#    because the mapping alignments are not always identical. Large gaps
#    might be inserted into features, causing the total feature length to grow.
#    Or regions covered by a feature might not be covered by alignments at all,
#    so the mapped feature is smaller than the original
#
#  --in_format <guess, hgvs, bed, gvf, gff, gtf, gff3, asnt, asnb, region>
#  --out_format
#  The input format of your annotation file, and the desired output format
#    Conversions to different formats might not preserve meta-data exactly.
#    Leave --out_format alone to get out the same format as in. 
#    Guess leaves the decision to our service.
#
#  --annot_out <local file name> 
#    The file this script will write the output annotation too, 
#    in the format from --out_format, or the same format as the --annotation file.
#  --report_out <local file name>
#    The file to write Remap's own report file to. It is a tab-delimited file,
#    and it includes base-by-base locations of where every segment moved through 
#    the alignments. It also explains why some segments do not map through.
#  --gbench_out <local file name>
#    The file to write a Gbench project to. This file will include the input 
#    and output annotation, and the alignments used.
#
#
#
my $fail_counter = 0;
my $MAX_FAILS = 5;

my $verbose = 0;
my $dev = 0;

my $help = 0;
my $mode = '';  # batches, asm-asm, asm-rsg, rsg-asm, alt-loci
my $annot_file_name = '';
my $from_acc = 'GCF_000001405.12';
my $dest_acc = 'GCF_000001405.13';

my $allow_dupes = 'on'; # on of off
my $merge = 'on'; # on or off
my $min_cov = 0.5; # intervals that remap less than this get dropped
my $max_exp = 2.0; # intervals that grow more than this get dropped

my $mapped_annot_file = '';
my $mapped_report_file = '';
my $mapped_gbench_file = '';

my $in_format = 'guess';
my $out_format = 'guess';

           GetOptions ( "help" => \$help,
                        "mode=s" => \$mode,                  # string
                        "annotation=s" => \$annot_file_name, # filename
                        "from=s" => \$from_acc,   # string, gencoll acc, or 'RefSeqGene'
                        "dest=s" => \$dest_acc,
                        "allowdupes=s" => \$allow_dupes,
                        "merge=s" => \$merge,
                        "mincov=s" => \$min_cov,
                        "maxexp=s" => \$max_exp,
                        "annot_out=s" => \$mapped_annot_file,
                        "report_out=s" => \$mapped_report_file,
                        "gbench_out=s" => \$mapped_gbench_file,
                        "in_format=s" => \$in_format,
                        "out_format=s" => \$out_format,
                        "verbose" => \$verbose,
                        "dev" => \$dev,
                        "jobid=s" => \$jobid,
                        "ctgtime=s" => \$ctgtime
                        ); 
                        
if(($help eq 1) or ($mode eq '')) {
    print( 
    " Usage: perl remap_api.pl
       
      Please limit your submissions to files of approximatly 250,000 rows, 
      and at most 4 simultaneous submissions.
      
      --mode <mode>
         batches, asm-asm, asm-rsg, rsg-asm, alt-loci.
       
      batches will print a list of all assembly-assembly pairs that alignments
      exist for.
      
      batches use no other argumentes, 
      other modes:
          asm-asm is for mapping between two assemblies
          asm-rsg is for mapping from an assembly to a gene data set (RefSeq or LRG)
          rsg-asm is for mapping from a gene data set to an assembly
          alt-loci is for mapping between a primary assembly and its related alt-loci

      required:
        --annotation <filename>
          The file containing your existing annotation (GFF3, GVF, BED, etc)
      
        --from <gencoll accession>
        --dest <gencoll accession>
          The gencoll accessions for the assemblies being mapped to and from
          Use 'RefSeqGene' when mapping to or from RefSeqGene
          Use 'LRG' when mapping to or from LRG
          Leave out --dest when use 'alt-loci' mode

      options:
        --allowdupes <'on' or 'off'>
          on (default) does all possible mappings,
          off excludes second pass mappings
      
        --merge <'on' or 'off'>
          on (default) will merge fragmented mappings back together
          off will leave fragmented results untouched
      
        --mincov <fraction 0.01 to 10.0, default 0.5 >
        --maxexp <fraction 0.01 to 10.0, default 2.0 >
          (r_cov = mapped_feature_length / original_feature_length
           if r_cov < min_cov then throw out feature
           if r_cov > max_exp then throw out feature )
          mapped_feature_length might not be identical to original_feature_length
          because the mapping alignments are not always identical. Large gaps
          might be inserted into features, causing the total feature length to grow.
          Or regions covered by a feature might not be covered by alignments at all,
          so the mapped feature is smaller than the original
      
        --in_format <guess (default), hgvs, bed, gvf, gff, gtf, gff3, asnt, asnb, region>
        --out_format
          The input format of your annotation file, and the desired output format
          Conversions to different formats might not preserve meta-data exactly.
          Leave --out_format alone to get out the same format as in.
          Guess leaves the decision to our service.  
        
        --annot_out <local file name> 
          The file this script will write the output annotation too, 
          in the format from --out_format, or the same format as the --annotation file.
        --report_out <local file name>
          The file to write Remap's own report file to. It is a tab-delimited file,
          and it includes base-by-base locations of where every segment moved through 
          the alignments. It also explains why some segments do not map through.
        --gbench_out <local file name>
          The file to write a Gbench project to. This file will include the input 
          and output annotation, and the alignments used.
        
          ");
    print("\n");

    exit(0);
}

if($verbose eq 1) {
    print "remap_api.pl"." : "."2016-12-22 : 16:44:12 "."\n";
} 
	
if($dev) {
    my $from = 'www';
    my $to = 'dev';
    $netloc =~ s/$from/$to/g;
    $URLbase = $protocol."://".$netloc."/".$remapdir;
}


my $userAgent = LWP::UserAgent->new(keep_alive=>1, ssl_opts => { verify_hostname => 1 });
# no credentials for Remap.
#$userAgent->credentials($netloc, $realm, $username, $password);
$userAgent->timeout(10);
$userAgent->env_proxy;


if($mode eq 'batches') {
    PrintBatches();
    exit(1);
}

elsif($mode eq 'asm-asm' or $mode eq 'asm-rsg' or $mode eq 'rsg-asm' or $mode eq 'alt-loci') {

    my $batchid = 0; # GetBatchId($from_acc, $dest_acc);
    
    if($mode eq 'asm-asm') {
        $batchid = GetAsmBatchId($from_acc, $dest_acc);
   
        if($batchid == 0) {
            print "Could not find an alignment batch for your assembly pair: " 
                 . $from_acc . " x " . $dest_acc . "\n";
            print "Please run \"--mode batches\" for a list of available assembly pairs.\n";
            exit(-1);
        }
    }

    my $main_result = 0;
    if($jobid eq "") {
        $main_result = MainSubmission($mode, $batchid, $from_acc, $dest_acc, $annot_file_name);
    } else {
        $forward_url = $URLbase.$RemapServiceCgi."?"."job_key=".$jobid."&ctg_time=".$ctgtime; 
    }

    if($main_result == 0) {
        print "Submitted: " . $jobid . "__" . $ctgtime . "\n";

        while(1) {
            sleep(3);
            my $recheck_result = RecheckSubmission($forward_url);
            if($recheck_result == 0) {
                print "Complete" . "\n";
                GetResultFile($annot_url, $mapped_annot_file);
                GetResultFile($report_url, $mapped_report_file);
                GetResultFile($gbench_url, $mapped_gbench_file);
                last;
            }
            else {
                print "Refreshing "."\n";
            }
        }
    }
    elsif($main_result == 1) {
        print "Some error happened. Your job was not submitted. \n";
        exit(-1);
    }
    else {
        print "Unexpected error: " . $main_result . "\n";
        exit(-1);
    }
}

else {
    print "Mode not recognized: " . $mode . "\n";
    exit(-1);
}

## end MAIN

# initial submission
# returns 0: Successful submission
# returns 1: Failed to send the submission was not accepted somehow
sub MainSubmission {

    my $mode = shift;
    my $batchid = shift;
    my $acc_one = shift;
    my $acc_two = shift;
    my $annot_file_name = shift;

    my $from_arg = 'source-assembly';
    my $dest_arg = 'target-assembly';

    if($mode eq 'rsg-asm' or $mode eq 'asm-rsg') {
        $from_arg = 'datamapfrom';
        $dest_arg = 'datamapto';
    }

    my $with_refseqgenes = 'on';
    my $without_refseqgenes = 'on';
    my $gbenchout = 'on';
    
    if ($mapped_gbench_file eq '') {
        $gbenchout = 'off';
    }

	# could provide an explicit version parameter. 
	# leave it out, and it returns 'latest'
	my $posturl = $URLbase . $RemapServiceCgi;


	my $request = (POST $posturl,
                    Content_Type => 'multipart/form-data',
                    'Keep-Alive' => 'timeout=60',
                    Content => [ mode => $mode , batch => $batchid ,
                                 $from_arg => $acc_one , 
                                 $dest_arg => $acc_two ,
                                 'allow_dupes' =>$allow_dupes , 'merge' => $merge , 
                                 'min_cov' => $min_cov , 'max_exp' => $max_exp ,
                                 'in_format' => $in_format , 'out_format' => $out_format ,
                                 'with_refseqgenes' => $with_refseqgenes ,
                                 'without_refseqgenes' => $without_refseqgenes ,
                                 'gbenchout' => $gbenchout , 
                                 'api' => 'true' ,
                                 annot => [$annot_file_name] ] ); 
    #if($verbose eq 1) {
    #    print $request->as_string . "\n";
    #} 
    my $response = $userAgent->simple_request($request);
    
    if($verbose eq 1) {
        print $posturl ."\n";
        print "MainSubmission"."\n";
        print $response->code . "\n";
        #print Dumper($response->headers()) . "\n";
        #print $response->content . "\n";
        print $response->as_string . "\n";
    }
    
    # HTTP reader NCBI-RCGI-RetryURL for retry url
    # parse XML out of response. jobid, jobstatus, ctgtime

    #<status>
    #<jobid>JSID_01_379_130.14.22.21_9001</jobid>
    #<job_status></job_status>
    #<ctgtime>1323790929</ctgtime>
    #<since>13 December 2011, 10:42:09</since>
    #<elapsed>0</elapsed>
    #<progress></progress>
    #<batch_id></batch_id>
    #</status>
	if($response->code == 403) {
        print "INIT RESP: 403"."\n";
        print "Probably too large a file. "."\n";
        return 1;
    } 
    elsif($response->code == 404) {
        print "INIT RESP: 404"."\n";
        print "probably too large a file."."\n";
        return 1;
    } 
    #if($response->is_success) {
		my $code = $response->code;
        $forward_url = $response->header("NCBI-RCGI-RetryURL");
        my $xml_content = $response->content;
        # this triggers an xml parse that includes connecting to DOCTYPEs and fetching them. 
        #   plz stop
        my $parser = XML::Parser->new( NoLWP => 1);
        my $xp = XML::XPath->new(parser => $parser, xml => $xml_content);
        $jobid = $xp->find("/status/jobid/text()");
        $ctgtime = $xp->find("/status/ctgtime/text()");
        my $jobstatus = $xp->find("/status/job_status/text()");
        
        #print "code: " . $code . "\n";
        #print "forward url: " . $forward_url . "\n";
        #print "xml body: " . $xml_content . "\n";
		#print "jobid: " . $jobid . "\n";
		#print "ctgtime: " . $ctgtime . "\n";
        #print "jobstatus: " . $jobstatus . "\n";

    if($code == 303) {
        return 0;
    } else {
        print "Error, MainSubmission: \n";
        print $response->as_string . "\n";
        return 1;
    }

}


# refresh caller function
# returns 1: Job is running, keep checking
# returns 0: Job is complete
# exit -2: Fatal error on server side
sub RecheckSubmission {

    my $recheck_url = shift;

    if($verbose eq 1) {
        print "recheck_url: ".$recheck_url . "\n";
    }

    my $request = (GET $recheck_url );

	my $response = $userAgent->request($request);

    if(not $response->is_success) {
        $fail_counter = $fail_counter + 1;
        if($fail_counter < $MAX_FAILS) {
            return 0;
        }

        print "Recheck of job ". $jobid."__".$ctgtime. " failed HTTP request.\n";
        print "This means Remap has crashed. It could be because the submitted data file was too large. Try splitting it up into smaller files? " ;
        
        if($verbose eq 1) {
            print "Recheck Fail"."\n";
            print "code: ".$response->code . "\n";
            print "header: ".$response->header . "\n";
            print "content: ".$response->content . "\n";
            print "as_string: ".$response->as_string . "\n";
        }
        exit(-2);
    }

    my $code = $response->code;
    my $curr_forward_url = $response->header("NCBI-RCGI-RetryURL");
    my $xml_content = $response->content;
   
    if($verbose eq 1) {
        print "recheck" . "\n"; 
        print "code: " . $code . "\n";
        print "header: ".$response->headers_as_string() . "\n";
        if(defined $curr_forward_url) {
            print "forward url: " . $curr_forward_url . "\n";
        }
        print "xml body: " . $xml_content . "\n";
        print "jobid: " . $jobid . "\n";
        print "ctgtime: " . $ctgtime . "\n";
    }

    if($code == 200 and defined($curr_forward_url) and $curr_forward_url ne '' ) {
        if($verbose eq 1) {
            print $response->as_string . "\n";
        }
        return 1;
    } else {
        # pull the links out of the response, the links
        # are the URLs to ncFetch the results
        my $parser = XML::Parser->new( NoLWP => 1);
        my $xp = XML::XPath->new(parser => $parser, xml => $xml_content);
        my $as = $xp->find("/HTML/BODY/div/span/a");
        
        foreach my $a ($as->get_nodelist) {
            my $href_url = "not_set";
            my $link_text = "not_set";
            foreach my $attrib ($a->getAttributeNodes()) {
                if( $attrib->getLocalName() eq "href" ) {
                    $href_url = $attrib->getNodeValue();
                }
            }
            foreach my $child ($a->getChildNodes()) {
                $link_text = $child->getData();
            }

            if($href_url ne "not_set") {
                if($link_text eq "Annot Link") {
                    $annot_url = $href_url;
                }
                elsif($link_text eq "Report Link") {
                    $report_url = $href_url;
                }
                elsif($link_text eq "GBench Project Link") {
                    $gbench_url = $href_url;
                }
                elsif($link_text eq "No features mapped to new locations.") {
                    print $link_text ."\n";
                    return 0;
                }
            }
            #print "FOUND\n\n", 
            #XML::XPath::XMLParser::as_string($a),
            #"\n\n";
        }
       
        # Failed to get a proper XML message, but don't fail abort,
        # just keep Rechecking
        if(!defined($annot_url) or $annot_url eq '') {
            # check for  <html><head><title>Empty Result :q
            #
            my $error_message = $xml_content;
            $error_message =~ s|<.+?>||g;
            #print $error_message ."\n";
            #print 'Failed'."\n";
          
            if( index($error_message, "Empty Result") != -1) {
                print 'An unknown error has occured on NCBI\'s side. Please contact us at info@ncbi.nlm.nih.gov '."\n"; 
                exit(-2);
            }

            if( index($error_message,"Error:") != -1 ) {
                print $error_message ."\n";
                exit(-2);
            }
            
            # unknown errors, keep trying
            return 1;
        }
    
        if($verbose eq 1) {
            print "Recheck Success Complete"."\n";
            print $response->code . "\n";
            print $response->as_string . "\n";
        }
     
        return 0;
    }

}


# Get Result File
sub GetResultFile {

    my $file_url = shift;
    my $save_file_name = shift;
    #print "* " . $file_url . " * \t * " . $save_file_name . " * \n";

    if($save_file_name eq '') {
        return 0;
    }

    if(!defined($file_url) or $file_url eq '') {
        print "Could not download file, empty URL given.";
        return 0;
    }
    
    my $local_fail_counter = 0;

    my $response;

    while(1) {
        my $request = (GET $file_url );
	    $response = $userAgent->request($request);

        if(not $response->is_success) {
            $local_fail_counter = $local_fail_counter + 1;
            if($local_fail_counter < $MAX_FAILS) {
                if($verbose eq 1) {
                    print "retrying download of ".$file_url."\n";
                }
                sleep(3);
                next;
            } else {
                print "download of '".$file_url."' failed.\n".$response->content ."\n";            
                return 0;
            }
        } else {
            last;
        }
    }

    my $code = $response->code;
    my $content = $response->content;
  
    #print $response->filename ;
    #print $response->as_string;

    if($save_file_name eq '') {
        $save_file_name = $response->filename;
    }

    unless(open SAVE, '>' . $save_file_name) {
        die "\nCannot create save file ".$save_file_name."\n";
    }
                       
    # Without this line, we may get a
    # 'wide characters in print' warning.
    #binmode(SAVE, ":utf8");
    print SAVE $response->content;
    close SAVE;
    print "\tSaving " . $save_file_name . "\n";
}


# function to use remapinfo.cgi, print info

sub PrintBatches {
    
	my $info_url = $URLbase . $RemapInfoCgi;
    
    if($verbose eq 1) {
        print "URL: ".$info_url."\n";
    }

    my $request = (GET $info_url );
	my $response = $userAgent->request($request);

    my $code = $response->code;
    my $json_content = $response->content;
 
    if($verbose eq 1) {
        print "Batches Verbose:"."\n";
        print "code: ".$response->code . "\n";
        print "header: ".$response->headers_as_string() . "\n";
        print "content: ".$response->content . "\n";
    }

    if(not $response->is_success) {
        print "Could not retrieve alignment batch list.\n";
        exit(-2);
    }

    $json_content = "{ \"RMP_Modes\": " .$json_content. " } ";

    #print $json_content;
    
    my $json = new JSON;
    my $fetched_json = $json->allow_nonref->utf8->relaxed->decode($json_content);

    my $RMPModes = $fetched_json->{"RMP_Modes"};
    my $RMPTaxlist = @$RMPModes[0]->{"RMP_Taxlist"};
    my $RMPAssemblylist = @$RMPModes[1]->{"RMP_Assemblylist"};
    my $RMPBatchlist = @$RMPModes[2]->{"RMP_Batchlist"};
    my $RMPRefseqgenelist = @$RMPModes[3]->{"RMP_Refseqgenelist"};
    my $RMPAltlocilist = @$RMPModes[4]->{"RMP_Altlocilist"};

    my $id_set = {  };
    
    my $tax_hr = {};
    my $assem_hr = {};
    foreach my $tax ( @$RMPTaxlist ) {
        my $curr = $tax;
        my $taxid = $curr->{"taxid"};
        $tax_hr->{ $taxid } = $curr;
    }
    foreach my $assem ( @$RMPAssemblylist ) {
        my $acc = $assem->{"accession"};
        $assem_hr->{ $acc } = $assem;
    }

    print "#". "batch_id" ."\t". 
        "query_species" ."\t".
        "query_name" ."\t".
        "query_ucsc" ."\t".
        "query_acc" ."\t".
        "target_species" ."\t".
        "target_name" ."\t".
        "target_ucsc" ."\t".
        "target_acc" ."\t".
        "alignment_date" ."\n" ;


    foreach my $batch ( @$RMPBatchlist ) {
        my $batch_id = $batch->{"batch_id"};
        my $query_acc = $batch->{query_acc};
        my $target_acc = $batch->{target_acc};
        my $datetime = $batch->{"create_date"};
        my $date = "";
        my $dtime= "";
        ($date, $dtime) = split(' ', $datetime, 2);


        if($id_set->{$batch_id}) {
            next;
        }
        $id_set->{$batch_id} = 1;

        my $query  = $assem_hr->{$query_acc};
        my $target = $assem_hr->{$target_acc};
        my $qtax = $tax_hr->{ $query->{taxid} };
        my $ttax = $tax_hr->{ $target->{taxid} };

        print $batch_id ."\t". 
            $qtax->{species} ."\t".
            $query->{name} ."\t".
            $query->{ucsc} ."\t".
            $query->{accession} ."\t".
            $ttax->{species} ."\t".
            $target->{name} ."\t".
            $target->{ucsc} ."\t".
            $target->{accession} ."\t".
            $date ."\n"
    }
  
    foreach my $assem ( @$RMPRefseqgenelist ) {
        my $batch_id = 0;
        
        my $tax = $tax_hr->{ $assem->{taxid} };

        print $batch_id ."\t". 
            $tax->{species} ."\t".
            $assem->{name} ."\t".
            $assem->{ucsc} ."\t".
            $assem->{accession} ."\t".
            $tax->{species} ."\t".
            "RefSeqGene" ."\t".
            "" ."\t".
            "RefSeqGene" ."\n"
    }

    foreach my $assem ( @$RMPAltlocilist ) {
        my $batch_id = 0;
        
        my $tax = $tax_hr->{ $assem->{taxid} };

        print $batch_id ."\t". 
            $tax->{species} ."\t".
            $assem->{name} ."\t".
            $assem->{ucsc} ."\t".
            $assem->{accession} ."\t".
            $tax->{species} ."\t".
            "AltLoci" ."\t".
            "" ."\t".
            "AltLoci" ."\n"
    }
}

# function to use remapinfo.cgi, and find a batchid
sub GetAsmBatchId {

    my $acc_one = shift;
    my $acc_two = shift;
    
	my $info_url = $URLbase . $RemapInfoCgi;
    my $request = (GET $info_url );

	my $response = $userAgent->request($request);

    my $code = $response->code;
    my $json_content = $response->content;
   
    $json_content = "{ \"RMP_Modes\": " .$json_content. " } ";

    #print $json_content;
    
    my $json = new JSON;
    my $fetched_json = $json->allow_nonref->utf8->relaxed->decode($json_content);

    my $RMPModes = $fetched_json->{"RMP_Modes"};
    my $RMPTaxlist = @$RMPModes[0]->{"RMP_Taxlist"};
    my $RMPAssemblylist = @$RMPModes[1]->{"RMP_Assemblylist"};
    my $RMPBatchlist = @$RMPModes[2]->{"RMP_Batchlist"};
    my $RMPRefseqgenelist = @$RMPModes[3]->{"RMP_Refseqgenelist"};
    my $RMPAltlocilist = @$RMPModes[4]->{"RMP_Altlocilist"};


    my $id_set = {  };
  
    foreach my $batch ( @$RMPBatchlist ) {
        my $batch_id = $batch->{"batch_id"};
        my $query_acc = $batch->{query_acc};
        my $target_acc = $batch->{target_acc};

        if($id_set->{$batch_id}) {
            next;
        }
        $id_set->{$batch_id} = 1;

        if( (($query_acc eq $acc_one) and
             ($target_acc eq $acc_two)) or
            (($query_acc eq $acc_two) and
             ($target_acc eq $acc_one)) ) {
            return $batch_id;
        }
    }

    return 0;
}



