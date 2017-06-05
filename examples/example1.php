<?php

/**
 * This script creates artificial (in silico) tetraploids that could be 
 * expected from autopolyploidy of the available diploid individuals.
 * 
 * In order to run the script, you need to place it, together with the 
 *  input file, into a folder inside document root of your http server. 
 * Then open your web browser and go to the address
 * http://{host}/{path/to/script}/example1.php
 * Where {host} is "localhost" when running on your local machine
 * Change {path/to/script} to actual path relative to the server's document root.
 * All information on Apache HTTP server can be found at http://httpd.apache.org/
 */

/**
 * Set here your input and output parameters
 */
//input file WITHOUT HEADERS AND LABELS
$inputFileName = 'AMAU_MSATdata_2x.txt';
//alleles in lines are separated by tabulator "\t", or space " "
$inputDelimiter = "\t";
//if data are microsatelites then true, else false
$isMSAT = true;
/* fill the number of alleles in each locus if applicable.
 * E.g. array(2, 4, 3) means that there are 3 loci, first locus contains
 * first 2 alleles, second locus contains next 4 alleles and the last 
 * locus contains next 3 alleles.
*/ 
$inputLoci = array(13, 12, 7, 6, 6, 3, 2, 2, 5, 6, 8);
//output file
$outputFileName = 'AMAU_MSATdata_2x';
//number of output files (sets of in silico tetraploids), e.g. 3
$outputFileCount = 3;
//number of in silico tetraploids in each file (randomly selected from all possible combinations), e.g. 30
$outputResultsCount = 30;

/**
 * Functions declarations
 */

/**
 * Reads input file, creates array of arrays of columns. 
 * Columns are separated by $delimiter.
 * Input file should be plain text data without headers.
 * @param string $filename input file
 * @param string $delimiter separator of 
 * @return array of lines
 */
function readFileToArray($filename, $delimiter) {
    if (empty($filename)) {
        die('Filename is empty');
    }
    $readF = fopen($filename, 'r') or die('Error opening file');
    $genArray = array();
    while ($line = fgets($readF)) {
        array_push($genArray, explode($delimiter, trim($line)));
    }
    fclose($readF);
    return $genArray;
}

/**
 * Creates combinations of additions of each line with every other line. Line
 * is not added to itself. Each pair is added once.
 * C = n! / ((n - 2)! * 2)
 * n: number of lines
 * @param array $lines array of arrays of alleles
 * @return array combinations
 */
function makeCombinations($lines) {
    if ($lines == NULL || empty($lines)) {
        die('Input array is null or empty');
    }
    $output = array();
    for ($j = 0; $j < count($lines) - 1; $j++) { //for each line
        $refline = $lines[$j]; //set current line as reference
        /*
         * for each line further in array, starting from the immediate next than the reference line
         */
        for ($i = $j + 1; $i < count($lines); $i++) {
            $plusline = $lines[$i]; //set line as the one to be added
            $res = addLines($refline, $plusline); //make addition
            array_push($output, $res); //add result to output
        }
    }
    return $output; //array of arrays
}

/**
 * Makes addition of two arrays of string. For the microsatelites, the array values
 * are either '1' or '0' - allele is present, allele is not present. The addition is 
 * performed on each position separately, formally it is a logical 'AND' operation.
 * If lineOne[position] == '1' and lineTwo[position] == '1' then result['position'] = '1'
 * else result[position] = 0
 * @param array $lineOne 
 * @param array $lineTwo
 * @return string
 */
function addLines($lineOne, $lineTwo) {
    if ($lineOne == NULL || empty($lineOne)) {
        die('First parameter is null or empty');
    }
    if ($lineTwo == NULL || empty($lineTwo)) {
        die('Second parameter is null or empty');
    }
    $result = array();
    for ($k = 0; $k < count($lineOne); $k++) {
        if (trim($lineOne[$k]) == "1" || trim($lineTwo[$k]) == "1") {
            $result[] = "1";
        } else {
            $result[] = "0";
        }
    }
    return $result;
}

/**
 * Check the input lines for following rule:
 * Nn any given group (locus), there must be at least one '1' (present allele) and 
 * the maximum count of '1's is given in group (locus) by threshold.
 * If a line does not hold this condition in any group, it is removed.
 * Group is a consecutive number of columns that belong together for any reason.
 * E.g. array(2, 4, 3) means that there are 3 loci, first locus contains
 * first 2 alleles, second locus contains next 4 alleles and the last locus contains next 3 alleles.
 * @param array $lines input array of lines to be tested
 * @param array $columnGroups specified groups (loci) of columns
 * @param integer $upperThreshold maximum number of '1's in one group
 * @return array lines that satisfy the condition on loci
 */
function checkCondition($lines, $columnGroups, $upperThreshold) {
    if ($lines == NULL || empty($lines)) {
        die('Input array si null or empty');
    }
    if ($columnGroups == NULL || empty($columnGroups)) {
        die('Loci array is null or empty');
    }
    $toremove = array();
    if (count($columnGroups) > 0) {
        foreach ($lines as $line) { //iterating each line
            $start = 0;
            foreach ($columnGroups as $num) { //for each locus
                $group = array_slice($line, $start, $num); //get alleles on a locus
                $occurences = array_count_values($group);
                //check for condition
                if (!in_array('1', $group) && $occurences['1'] > $upperThreshold) {
                    array_push($toremove, $line); //if holds, add line to removed
                    break;
                }
                $start += $num;
            }
        }
    }
    //substract removed lines from input lines
    $final = array_diff($lines, $toremove);
    return $final;
}

/**
 * Creates number of specified files. For each file, random lines are chosen such that they
 * do not repeat. Number of chosen lines is given by parameter.
 * @param array $lines input array of lines to choose from
 * @param integer $filesCount number of files to be generated
 * @param integer $resultsCount number of results to be picked
 * @param string $outfile prefix of output file. The resulting 
 * files are suffixed with '_results_i.txt' where i is increment of number of files
 */
function pickRandom($lines, $filesCount, $resultsCount, $outfile) {
    if ($lines == NULL || empty($lines)) {
        die('Input array is null or empty');
    }
    if (!is_int($filesCount) || $filesCount <= 0) {
        die('Files count must be integer greater than 0');
    }
    for ($i = 0; $i < $filesCount; $i++) { //repeat for each file we want to create
        $f = $i + 1;
        $writeF = fopen($outfile . "_results_$f.txt", 'w') or die('Error opening file for writing');
        $chosen = array(); //array of chosen lines, only numbers of lines
        for ($j = 0; $j < $resultsCount; $j++) { //reapeat as many times as is resultsCount
            $rand = rand(0, count($lines) - 1); //generate random number
            while (in_array($rand, $chosen)) { //if the line number has been chosen before
                $rand = rand(0, count($lines) - 1); //generate new random number until it is a new one
            }
            array_push($chosen, $rand); //add line number to chosen
            $out_str = implode("\t", $lines[$rand]); //get actual line and join characters with delimiter
            fwrite($writeF, $out_str . "\n"); //write to file
        }
        fclose($writeF);
    }
}

/**
 * Execution part
 */
//read file
$lines = readFileToArray($inputFileName, $inputDelimiter);
echo "Loaded " . count($lines) . " lines<br />\n";
//create full set
$combinations = makeCombinations($lines);
echo "Created " . count($combinations) . " combinations<br />\n";
//check loci conditions
if ($isMSAT) {
    $linesClean = checkCondition($combinations, $inputLoci, 4);
    echo count($linesClean) . " lines remaining after checking condition<br />\n";
} else {
    $linesClean = $combinations;
}
//choose random results
echo "Creating result files<br />\n";
pickRandom($linesClean, $outputFileCount, $outputResultsCount, $outputFileName);
echo "Finished";
