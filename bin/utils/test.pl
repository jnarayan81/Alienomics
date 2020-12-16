sub toString
    {
       my $self = shift;
       my $thishash = $self->{_listofkeys};
       my $memfile;
       open (FORMATTED_OUT, '>',\$memfile) || warn "unable to open memory file!";

format FORMATTED_OUT =
        IMEI              BrandName            Model             OS                  Registered To
        ==========================================================================================================
        @<<<<<<<<<<<<<<<< @<<<<<<<<<<<<<<<<<<< @<<<<<<<<<<<<<<<< @<<<<<<<<<<<<<<<<<< @<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        $thishash->{IMEI}, $thishash->{PhoneBrandName}, $thishash->{PhoneModel}, $thishash->{PhoneOS}, $thishash->{RegisteredTo}
.

    write FORMATTED_OUT ;
    close (FORMATTED_OUT);

print $memfile;
    return  $memfile;
    }



toString($ARGV[0]);
