cmd=paste(
  "my $dir=shift @ARGV;",
  "open DATA,\"$dir/temp_files/read_nonsplit.txt\" or die;",
  "my %chrs;",
  "<DATA>;",
  "while(my $str=<DATA>){",
  "    my @str=split(/\t/,$str);",
  "    if(defined($chrs{$str[3]})){",
  "        print {$chrs{$str[3]}} \"$str[0]\t$str[4]\t$str[5]\";",
  "    }",
  "    else{",
  "      open(my $fileout,\">\",\"${dir}/temp_files/read_nonsplit_$str[3].txt\");",
  "      print $fileout \"$str[0]\t$str[4]\t$str[5]\";",
  "      $chrs{$str[3]}=$fileout;",
  "    }",
  "}",
  "close DATA;",

  "foreach (keys %chrs){",
  " close($chrs{$_});",
  "}",
  sep="\n"
)

