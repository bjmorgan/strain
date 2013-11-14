#!/usr/local/bin/ruby

# fstrain, Benjamin J. Morgan, 2013.

# Calculates local strain fields from VASP / QUAIM ionic configurations
# This routine is based on the algorithm proposed by 
# M. L. Falk  and J. S. Langer; Phys. Rev. E 57, 7192 (1998),
# and borrows heavily from the MATLAB implementation by P. Schall and D. Bonner
# http://www.science.uva.nl/research/wzi/scm/strainfield.m
# used in, for example, Rahmani et al. Sci. Rep. 2, 1064 (2013).

# Three output files are produced:
# epsilon.dat => 

require 'matrix'

class Ion

	attr_reader :r, :number, :epsilon

	@@ion_number = 0

  def initialize( r )
  	@@ion_number += 1
  	@number = @@ion_number
    @r = r
    @neighbour_vectors = []
    @neighbours = []
    @epsilon = nil
  end

  def find_neighbours( ions, cell, cutoff )
  	@neighbours = []
  	@neighbour_vectors = []
  	ions.reject{ |i| i == self }.each do |ion|
  		dr = cell.dr( ion.r, @r )
  		@neighbour_vectors << dr if Vector[ *dr ].magnitude < cutoff
  	end
  end

  def strain_vector( reference_vectors, cell )
  	calculate_strain_epsilon( reference_vectors, cell ) if @epsilon.nil?
  	e = Matrix[ *@epsilon ].eigen
  	eigenvectors = e.v.column_vectors
  	eigenvalues = (0..2).collect{ |i| e.d[i,i] }
  	max_eigenvector = eigenvectors.zip( eigenvalues ).max_by{ |vector, value| value }
  	( max_eigenvector[1] * max_eigenvector[0] ).to_a
  end

  def calculate_strain_epsilon( reference_vectors, cell )
  	zero_strain_vectors = @neighbour_vectors.map{ |v| closest_reference_vector( v, reference_vectors, cell ) }
  	x = falk_x( zero_strain_vectors, @neighbour_vectors )
  	y = falk_y( zero_strain_vectors )
  	deformation_matrix = Matrix[ *affine_deformation( x, y ) ]
  	symmetric_component = ( 0.5 * ( deformation_matrix + deformation_matrix.transpose ) ).to_a
  	@epsilon = symmetric_component
  end

  def d_squared( reference_vectors, cell )
  	calculate_strain_epsilon( reference_vectors, cell ) if @epsilon.nil?
  	zero_strain_vectors = @neighbour_vectors.map{ |v| closest_reference_vector( v, reference_vectors, cell ) }
  	zero_strain_vectors.zip( @neighbour_vectors ).inject(0.0) do |sum_n, (v1,v2)| # loop over index n
  		sum_n + v1.each_with_index.inject(0.0) do |sum_i, (r_i, i)| 								# loop over index i
  			sum_i + ( r_i - v2.each_with_index.inject(0.0) do |sum_j, (r_j, j)| 			# loop over index j
  				sum_j + ( @epsilon[i][j] + delta(i,j) ) * r_j 
	  		end )**2
			end
  	end
  end

  def delta(i,j)
  	i == j ? 1.0 : 0.0
  end

  def falk_x( zero_strain_vectors, neighbour_vectors )
  	x = [ [ 0.0, 0.0, 0.0], [ 0.0, 0.0, 0.0 ], [ 0.0, 0.0, 0.0 ] ]
  	x.each_with_index do |row,i|
  		row.each_index do |j|
  				zero_strain_vectors.zip( neighbour_vectors ).each{ |v1,v2| x[i][j] += v1[i] * v2[j] }
  		end
  	end
  end

  def falk_y( zero_strain_vectors )
  	y = [ [ 0.0, 0.0, 0.0], [ 0.0, 0.0, 0.0 ], [ 0.0, 0.0, 0.0 ] ]
  	y.each_with_index do |row,i|
  		row.each_index do |j|
  				zero_strain_vectors.each{ |v| y[i][j] += v[i] * v[j] }
  		end
  	end
  end

  def affine_deformation( x, y )
  	y_inv = Matrix[ *y ].inv.to_a
  	deform = [ [ 0.0, 0.0, 0.0], [ 0.0, 0.0, 0.0 ], [ 0.0, 0.0, 0.0 ] ]
  	deform.each_with_index do |row,i|
  		row.each_index do |j|
  			3.times do |k|
  				deform[i][j] += x[i][k] * y_inv[j][k]
  			end
  		end
  	end
  	deform = ( Matrix[ *deform ] - Matrix.identity(3) ).to_a
  	deform
  end

  def closest_reference_vector( neighbour_vector, reference_vectors, cell ) 
  	reference_vectors.min_by { |vector| Vector[ *cell.dr( vector, neighbour_vector ) ].magnitude }
  end

end

class Poscar

	attr_reader :lattice, :coordinates

	def initialize( title, scaling, lattice, species, coordinate_type, coordinates )
		@title = title
		@scaling = scaling
		@lattice = Cell.new( lattice )
		@species = species
		@coordinate_type = coordinate_type
		@coordinates = coordinates
	end

	def self.new_from_file( filename, scaled, coordinate_type )
		abort( "VASP file #{filename} not found" ) unless File.file?( filename )
		poscar = File.new( filename, 'r' ).readlines
		title = poscar.shift.strip
		scaling = poscar.shift.strip.to_f
		lattice = poscar.shift(3).map{ |line| line.split.map{ |s| s.to_f } }
		lattice.map!{ |line| line.map!{ |f| f * scaling } } if scaled
		this_cell = Cell.new( lattice )
		species_names = poscar.shift.split
		species_numbers = poscar.shift.split.map{ |s| s.to_i }
		coordinates_type = poscar.shift
		total_ions = species_numbers.inject(0) { |result, element| result + element }
		coordinates = poscar.shift( total_ions ).map { |line| line.split.map { |s| s.to_f } }
		coordinates.map!{ |r| r.map! { |f| f * scaling } } if ( coordinates_type.match( /^[cC]/ ) and scaled )
		case coordinate_type
		when :cartesian
			coordinates.map!{ |r| this_cell.direct_to_cartesian_coordinates( r ) } if coordinates_type.match( /^[dD]/ )
		when :direct
			coordinates.map!{ |r| this_cell.cartesian_to_direct_coordinates( r ) } if coordinates_type.match( /^[cC]/ )
		else
			raise ( "\`#{coordinate_type}\' is not a valid coordinate type" )
		end
		Poscar.new( title, 1.0, lattice, Hash[ species_names.zip( species_numbers ) ], :cartesian, coordinates )
	end 
end

class Cell

  def initialize( l )
  	@lattice = l
  	@cell_lengths = @lattice.map{ |v| Vector[ *v ].magnitude }
  end

	def direct_to_cartesian_coordinates( r )
		( Matrix[ *@lattice ] * Matrix.column_vector( r ) ).column(0).to_a
	end

	def cartesian_to_direct_coordinates( r )
		( Matrix[ *@lattice ].inv * Matrix.column_vector( r ) ).column(0).to_a
	end

	def cartesian_minimum_image( vector )
	 	return_vector = cartesian_to_direct_coordinates( vector ).to_a.map do |r| 
	 		r -= 1.0 if r > 0.5
	 		r += 1.0 if r < -0.5 
	 		r
	 	end 
	 	return direct_to_cartesian_coordinates( return_vector )
  end

  def dr( r1, r2 )
  	cartesian_minimum_image( ( Vector[ *r1 ] - Vector[ *r2 ] ).to_a )
  end

end

class Calculation_Type

	attr_reader :default_positions_filename, :default_cell_filename

	def initialize( string, re, default_positions_filename, default_cell_filename )
		@string = string
		@re = re
		@default_positions_filename = default_positions_filename
		@default_cell_filename = default_cell_filename
	end

	def match( string )
		return @re.match( string )
	end

end

calculation_types = { :vasp => Calculation_Type.new( 'vasp', /^[vV]/, 'POSCAR', nil ),
											:quaim => Calculation_Type.new( 'quaim', /^[qQ]/, 'poscart.out', 'cellbox.out' ) }

require 'optparse'

options = {}

executable_name = File.basename( $PROGRAM_NAME )

optparse = OptionParser.new do |opts|
  opts.banner = "Usage: #{executable_name} [options] data_filenames"

# Set input type (currently VASP or QUAIM)
 	options[:input_type] = nil
 	opts.on( '-i', '--input TYPE', 'Input type (vasp|quaim)' ) do |type|
 		this_type = calculation_types.select{ |k,v| v.match( type ) }.keys[0]
 		if this_type.nil?
 		 	abort( "\"#{type}\" is not a valid input type\n" )
 		else
	 		options[:input_type] = this_type
	 	end
  end

# Read in positions filename
  options[ :positions_filename ] = nil
  opts.on( '-p', '--positions_file FILE', 'Positions file' ) do |file|
  	abort( "Positions file #{file} not found" ) unless File.file?( file )
  	options[ :positions_filename ] = file
  end

# Read in cell box filename
  options[ :box_filename ] = nil
  opts.on( '-b', '--box_file FILE', 'Box lattice file' ) do |file|
  	abort( "Box file #{file} not found" ) unless File.file?( file )
  	options[ :box_filename ] = file
  end

# Read in neighbour vectors filename
  options[ :vectors_filename ] = nil
  opts.on( '-v', '--vectors_file FILE', 'Neighbour vector file' ) do |file|
  	abort( "Neighbour vector file #{file} not found" ) unless File.file?( file )
  	options[ :vectors_file ] = file
  end

  options[ :cutoff ] = 2.5
  opts.on( '-c', '--cutoff CUTOFF', Float, "Neighbour cutoff distance (default #{options[ :cutoff ]})" ) do |cutoff|
  	options[ :cutoff ] = cutoff
  end

  options[ :scaling ] = 1.0
  opts.on( '-s', '--vector_scaling SCALE', Float, "Scale neighbour vectors by SCALE (default #{options[ :scaling ]})" ) do |scale|
  	options[ :scaling ] = scale
  end

  options[ :vector_scaling ] = 1.0
  opts.on( '-a', '--vector_scaling SCALE', Float, "Scale vector length in output by SCALE (default #{options[ :vector_scaling ]}" ) do |scale|
  	options[ :vector_scaling ] = scale
  end

end

def read_vasp_data_from_file( filename )
	poscar = Poscar.new_from_file( filename, true, :cartesian )
	return poscar.lattice, poscar.coordinates
end

def read_quaim_data_from_files( positions_filename, box_filename )
	abort( "quaim not implemented yet")
end

def read_vector_data_from_file( filename )
	File.new( filename, 'r' ).readlines.map{ |line| line.split.map{ |s| s.to_f } }.reject{ |line| line.empty? }
end

optparse.order!

# Use filename defaults if none have been provided by the user
positions_filename = options[:positions_filename] ||= calculation_types[ options[:input_type] ].default_positions_filename
box_filename       = options[:box_filename]       ||= calculation_types[ options[:input_type] ].default_cell_filename
vectors_filename   = options[:vectors_filename]   ||= 'neighbour_vectors.xyz'

case options[:input_type]
when :vasp
	cell, coordinates = read_vasp_data_from_file( positions_filename )
when :quaim
	cell, coordinates = read_quaim_data_from_files( positions_filename, box_filename )
else
	abort( "error parsing input type" )
end

ions = coordinates.map{ |c| Ion.new( c ) }

neighbour_vectors = read_vector_data_from_file( vectors_filename ).map{ |v| v.map{ |f| f * options[:scaling] } }

def output_format( n )
	return Array.new( n, '%.3e' ).join(' ')
end

epsilon_out   = File.new( 'epsilon.dat', 'w' )
vector_out    = File.new( 'eigenvector.dat', 'w' )
d_squared_out = File.new( 'dsquared.dat', 'w' )

epsilon_out.puts   "\# n  eps_xx     eps_xy     eps_xz     eps_yx     eps_yy     eps_yz     eps_zx     eps_zy     eps_zz"
vector_out.puts    "\# n  x          y          z          dx         dy         dz"
d_squared_out.puts "\# n  D^2"

ions.each_with_index do |ion,i| 
	ion.find_neighbours( ions, cell, options[:cutoff] )
	ion.calculate_strain_epsilon( neighbour_vectors, cell )
	strain = ion.strain_vector( neighbour_vectors, cell ).map{ |r|  r * options[:vector_scaling] }
	d2 = ion.d_squared( neighbour_vectors, cell )

	epsilon_out.puts   "  #{i+1}  #{output_format(9) % ion.epsilon.flatten}"
	vector_out.puts    "  #{output_format(3) % ion.r}  #{output_format(3) % strain }"
	d_squared_out.puts "  #{i+1}  #{output_format(1) % d2}}"
end
