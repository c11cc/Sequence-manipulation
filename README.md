# Sequence-manipulation

`gc_data.plx`<br>
This is a script for calculating the GC content of settled length around target site at a box length continuously (from the first base to the last base of the settled region with the box length as counting unit);<br>
* desired site infomation file contains chromosome and site infomation in the first and second column,no headings; <br>
* desired sequence file is in .fa format, with chrmosome id after >;<br>
* desired region should be smaller than the length of the sequence;<br>
* GC box length is the length of the region for calculating GC content;<br>
* base 'N' is excluded from calculation.<br>

`seq_extract.plx`<br>
This is a script for fetching the sequence of settled length around target site from fasta file;<br>
* desired site infomation file contains identifier and site infomation in the first and second column,no headings;
* desired sequence file is in .fa format, with identifier after >;<br>
* desired region should be smaller than the length of the sequences.
