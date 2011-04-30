#!/usr/bin/env ruby

# Run migrate-n a certain number of times and average runtime performance.
# Syntax: run.rb [--mp] [times]

require 'fileutils'
require 'benchmark'

EXAMPLES_DIR = File.expand_path(File.join(File.dirname(__FILE__), '..', 'example'))
PARMFILE = 'parmfile.testbayes'
SP_RUN = "../src/migrate-n #{PARMFILE} -nomenu"
MP_RUN = "mpirun -c 4 ../src/migrate-n-mpi #{PARMFILE} -nomenu"

def run
  times = 4
  command = SP_RUN
  if ARGV.first == "--mp"
    command = MP_RUN
    ARGV.shift
  end
  if ARGV.first != nil
    times = ARGV.first.to_i
  end
  
  puts "Running #{command} #{times} times"
  
  FileUtils.cd(EXAMPLES_DIR)
  results = (1..times).map do
    print "Running ... "
    time = Benchmark.realtime do
      system "#{command} &> benchmark.log"
      if $? != 0
        raise "Last run gave return code of #{$?}, exiting."
      end
    end
    puts "took #{time} seconds"
    time
  end
  
  puts results
end

run