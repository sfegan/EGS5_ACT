#!/usr/bin/perl -w

#
# Program: code_gen_array_ini_cpp.pl
# Author:  Stephen Fegan
#          Physics and Astronomy Department, UCLA
# E-mail:  sfegan@astro.ucla.edu
# Date:    2004/12/06
#
# The code generates reads in the file ARRAY_INI_TEMPLATE and generates
# the VSOArrayParameters.cpp code file.
#
# These two code generators read in the entries in an "array.ini" file
# (or the ARRAY_INI_TEMPLATE) and produce the C++ code required to
# read and write an array configuration to/from the "array.ini" file
# and transfer these values to/from the database. Initially the C++
# code was written by hand but it became very difficult to add values
# to or modify the structure of the "array.ini" file because the C++
# code would have to be changed by hand at the same time. This was
# very inconvenient and left the code prone to errors, so this code
# generator was developed. Hopefully this system will be relatively
# easy to maintain, although it does require the knowledge of PERL.
#
# The code proceeds as follows.
#
# (1) The file is read and the comments/members/units/types stored
#
# (2) The .hpp and .cpp files are generated using the stored data
#
# The ARRAY_INI_TEMPLATE should have entries as follows:
#
#   Comment comment comment comment comment comment comment comment
#   comment comment comment comment comment comment comment comment
# @ NameOfVariable [Units] <Variable Type>
#   Value
#
# For example:
#
#   Vector to intersection of rotation axes from the center of
#   reflector in reflector reference frame (+y points along optical axis,
#   +z is up, and +x points East when telescope is in Home position,
#   the origin is in the center of the reflector).
# @ ScopeTranslationX [cm] <double>
#   10.0
#

use strict;

my @Comments;
my @Members;
my @Types;
my @Units;
my @Values;

my @comment_bits;
my $ampersand=0;

foreach (<ARGV>)
  {
    chomp;
    s/\s+$//;

    if($ampersand==1)
      {
	$ampersand=0;
	s/^\s+//;
	s/\s+$//;
	s/^"//;
	s/"$//;
	push @Values,$_;
	next;
      }

    next unless $_;

    push(@comment_bits,$_);

    if(/^\@/)
      {
	my ( $amp, $var, $rest ) = split /\s+/,$_,3;
	my $unit="";
	my $type="double";
	
	$unit = $1 if($rest =~ /\[(.*)\]/);
	$type = $1 if($rest =~ /<(.*)>/);

	push @Comments,join("\\n",@comment_bits);
	push @Members,$var;
	push @Types,$type;
	push @Units,$unit;

	undef @comment_bits;

	$ampersand = 1;
      }
  }

my @comments;
my @members;
my @types;
my @units;
my @values;

my $comment;
my $member;
my $first;
my $line;

# =============================================================================
# =============================================================================
# =============================================================================

print << 'END';
//-*-mode:c++; mode:font-lock;-*-

// AUTOMATICALLY GENERATED by code_gen_array_ini_cpp.pl DO NOT EDIT!!

#include <fstream>
#include <iostream>
#include <sstream>
#include <memory>
#include <vector>
#include <limits>
#include <cassert>

#include <VSDataConverter.hpp>

#include <VSLineTokenizer.hpp>

#include "VSOArrayParameters.hpp"

using namespace Physics;
using namespace VERITAS;

const std::string VSOArrayParameters::scCollection("ArrayINI");

#define inf std::numeric_limits<double>::infinity();
END

# =============================================================================
# =============================================================================
# =============================================================================

@members=@Members;
@types=@Types;
@values=@Values;
while($member = shift @members)
  {
    my $type = shift @types;
    my $val = shift @values;

    printf("%-10s VSOArrayParameters::sCanonical%-30s = VSReturnedTypeDatumConverter<%s>::fromString(\"%s\");\n",$type,$member,$type,$val);
  }

print << 'END';

VSOArrayParameters::VSOArrayParameters(const VSOArrayParameters& o):
END

# =============================================================================
# =============================================================================
# =============================================================================

@comments=@Comments;
@members=@Members;
@types=@Types;
@units=@Units;

$first=1;
$line=" ";
while($member = shift @members)
  {
    $line .= "," unless $first;
    if(length($line)+3+length($member) > 78)
      {
	print $line,"\n";
	$line = " ";
      }

    $line .= " ".$member."(o.".$member.")";
    $first=0;
  }
print $line,"\n";

# =============================================================================
# =============================================================================
# =============================================================================

print << 'END';
{
  // nothing to see here
}

VSOArrayParameters::VSOArrayParameters(bool use_canonical_values):
END

# =============================================================================
# =============================================================================
# =============================================================================

@comments=@Comments;
@members=@Members;
@types=@Types;
@units=@Units;

$first=1;
$line=" ";
while($member = shift @members)
  {
    $line .= "," unless $first;
    if(length($line)+3+length($member) > 78)
      {
	print $line,"\n";
	$line = " ";
      }

    $line .= " ".$member."()";
    $first=0;
  }
print $line,"\n";

# =============================================================================
# =============================================================================
# =============================================================================

print << 'END';
{
  if(use_canonical_values)
    {
END

# =============================================================================
# =============================================================================
# =============================================================================

@comments=@Comments;
@members=@Members;
@types=@Types;
@units=@Units;
@values=@Values;

$first=1;
while($member = shift @members)
  {
    print "      ",sprintf("%-25s",$member)," = sCanonical",$member,";\n";
  }

# =============================================================================
# =============================================================================
# =============================================================================

print << 'END';
    }
}

VSOArrayParameters::~VSOArrayParameters()
{
  // nothing to see here
}

void VSOArrayParameters::dump(std::ostream& stream)
{
END

# =============================================================================
# =============================================================================
# =============================================================================

@comments=@Comments;
@members=@Members;
@types=@Types;
@units=@Units;
@values=@Values;

while($member = shift @members)
  {
    print "  stream << \"",sprintf("%-25s",$member)," = \" << ",$member," << std::endl;\n";
  }

# =============================================================================
# =============================================================================
# =============================================================================

print << 'END';
}

void VSOArrayParameters::reset(bool use_canonical_values)
{
  if(use_canonical_values)
    {
END

# =============================================================================
# =============================================================================
# =============================================================================

@members=@Members;
while($member = shift @members)
  {
    printf ("      %-25s = sCanonical%s;\n",$member,$member);
  }
print "    }\n  else\n    {\n";
@members=@Members;
@types=@Types;
while($member = shift @members)
  {
    printf ("      %-25s = %s();\n",$member,shift @types);
  }

# =============================================================================
# =============================================================================
# =============================================================================

print << 'END';
    }
}

static const char* array_ini_comment[] = {
END

# =============================================================================
# =============================================================================
# =============================================================================

@comments=@Comments;
@members=@Members;
@types=@Types;
@units=@Units;

$first=1;
while($comment = shift @comments)
  {
    print ",\n" unless $first;
    print "  \"",$comment,"\"";
    $first = 0;
  }
print "\n";

# =============================================================================
# =============================================================================
# =============================================================================

print << 'END';
};

#define NARRAYINIITEMS (sizeof(array_ini_comment)/sizeof(*array_ini_comment))

void VSOArrayParameters::zeroCanonicalValues()
{
END

# =============================================================================
# =============================================================================
# =============================================================================

@members=@Members;
@types=@Types;
while($member = shift @members)
  {
    printf("  sCanonical%-25s = %s();\n",$member,shift @types);
  }

# =============================================================================
# =============================================================================
# =============================================================================

print << 'END';
}

void VSOArrayParameters::resetCanonicalValues()
{
END

# =============================================================================
# =============================================================================
# =============================================================================

@members=@Members;
@types=@Types;
@values=@Values;
while($member = shift @members)
  {
    printf("  sCanonical%-25s = VSReturnedTypeDatumConverter<%s>::fromString(\"%s\");\n",$member,shift @types,shift @values);
  }

# =============================================================================
# =============================================================================
# =============================================================================

print << 'END';
}

void VSOArrayParameters::setCanonicalValuesFromArrayParameters(const VSOArrayParameters& o)
{
END

# =============================================================================
# =============================================================================
# =============================================================================

@comments=@Comments;
@members=@Members;
@types=@Types;
@units=@Units;

while($member = shift @members)
  {
    printf("  sCanonical%-25s = o.%s;\n",$member,$member);
  }

# =============================================================================
# =============================================================================
# =============================================================================

print << 'END';
}

bool VSOArrayParameters::
readFromArrayINIFile(const std::string& filename)
{
  std::ifstream stream(filename.c_str());
  if(stream)return readFromArrayINIFile(stream);
  else return false;
}

bool VSOArrayParameters::
readFromArrayINIFile(std::istream& stream)
{
  VSLineTokenizer tokenizer;

  std::vector<VSToken> parameters;
  while(!stream.eof())
    {
      VSTokenList tokens;
      tokenizer.tokenize(stream,tokens);

      if((tokens.size()!=0)&&(tokens[0].string().size()>0)&&
         (tokens[0].string().at(0)=='@'))
	{
          tokenizer.tokenize(stream,tokens);
          if(tokens.size() != 0)parameters.push_back(tokens[0]);
	}
    }

  if(parameters.size() != NARRAYINIITEMS)return false;

  unsigned p=0;
  bool good = true;

END

# =============================================================================
# =============================================================================
# =============================================================================

@comments=@Comments;
@members=@Members;
@types=@Types;
@units=@Units;

while($member = shift @members)
  {
    print "  good &= parameters[p++].convertTo(",$member,");\n";
  }

# =============================================================================
# =============================================================================
# =============================================================================

print << 'END';

  return((p==NARRAYINIITEMS)&&(good));
}

void VSOArrayParameters::
writeToArrayINIFile(const std::string& filename) const
{
  std::ofstream stream(filename.c_str());
  writeToArrayINIFile(stream);
}

void VSOArrayParameters::
writeToArrayINIFile(std::ostream& stream) const
{
  std::vector<std::string> parameters;

END

# =============================================================================
# =============================================================================
# =============================================================================

@comments=@Comments;
@members=@Members;
@types=@Types;
@units=@Units;

while($member = shift @members)
  {
    print "  parameters.push_back(VSDataConverter::toString(",$member,"));\n";
  }

# =============================================================================
# =============================================================================
# =============================================================================

print << 'END';

  assert(parameters.size() == NARRAYINIITEMS);

  for(unsigned i=0;i<parameters.size();i++)
    {
      if(i)stream << std::endl;
      stream << array_ini_comment[i] << std::endl
             << "  " << parameters[i] << std::endl;
    }
}

END

# =============================================================================
# =============================================================================
# =============================================================================

print << 'END';

#ifdef TEST_MAIN
int main(int argc, char** argv)
{
  VSOArrayParameters param;
  if(param.readFromArrayINIFile("array.ini"))
    param.writeToArrayINIFile(std::cout);
}
#endif

#ifdef TEST_MAIN2
int main(int argc, char** argv)
{
  Options options(argc,argv);
  VSOArrayParameters::zeroCanonicalValues();
  VSOArrayParameters::setCanonicalValuesFromOptions(options);
  VSOArrayParameters param(true);
  param.dump(std::cout);
  for(unsigned i=0;i<argc;i++){ if(i!=0)std::cout << ' '; std::cout << argv[i]; }
  std::cout << std::endl;
}
#endif
END
