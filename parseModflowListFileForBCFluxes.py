"""
    The intended use of this program is to support development of
    the North Florida-Southeast Georgia Groundwater Flow Model.
    I've done some error checking, but can't guarantee that I found
    checked everything, so please at your own risk and report any
    any errors to Trey Grubbs (jwg@srwmd.org).
"""
import time
import shelve

class GagedReaches:
    """ Class for working with gaged-reach data. """
    def __init__(self, reach_definition_file_name):
        self.define_gaged_reaches_by_bc_reach_indices(reach_definition_file_name)

    def define_gaged_reaches_by_bc_reach_indices(self, reach_definition_file_name):
        """ Read combinations of gaged-reach identifier and boundary-condition
           'reach indeces' that define each gaged reach. The boundary-condition
            reach indices correspond to to from MODFLOW listing file) data that define each gaged reach.
        
            The parameter, reach_definition_file_names, is a string containing
            the name of the file defining the reach associations for River Package
            Drain Package records.
            This method expects comma-delimited records to have the following fields:
            
            gaged_reach_id:             a string identifying the reach
            bc_type:                    a string indicating the reach type. Currently implmented
                                        bc types include 'drn' (for MODFLOW Drain Package),
                                        'ghb' (for MODFLOW General Head Boundary Package)
                                        'riv' for (MODFLOW River Package)
            bc_reach_id:                an integer that corresponds to the reach number
                                        in the table of simulated River Package leakage
                                        values reported in the MODFLOW listing file
        """

        self.gaged_reach_list = []
        self.gaged_reach_bc_reach_types = {}
        self.gaged_reach_bc_reach_list = {}

        f = open(reach_definition_file_name, 'r')

        for line in f.readlines():
            line_list = line.rstrip().split(',')
            # convert bc reach identifiers to integers
            line_list[1] = int(line_list[1])
            gaged_reach_id, bc_reach_id, bc_type = tuple(line_list)
            if self.gaged_reach_bc_reach_list.has_key((gaged_reach_id, bc_type)):
                if (bc_reach_id in self.gaged_reach_bc_reach_list[(gaged_reach_id, bc_type)]) == False:
                    self.gaged_reach_bc_reach_list[(gaged_reach_id, bc_type)].append(bc_reach_id)
            else:
                if gaged_reach_id not in self.gaged_reach_list:
                    # first input record for this gaged reach
                    self.gaged_reach_list.append(gaged_reach_id)
                    self.gaged_reach_bc_reach_types[(gaged_reach_id)] = [bc_type]
                else:
                    # first occurrence of this bc_type for this gaged_reach_id
                    self.gaged_reach_bc_reach_types[(gaged_reach_id, 'bc_type_list')].append(bc_type)
                self.gaged_reach_bc_reach_list[(gaged_reach_id, bc_type)] = [bc_reach_id]

        f.close()

    def assign_observations_to_reaches(self, reach_observations_file_name):
        """ Read a file containing observed values for each reach.
        
            This file is specified in the input parameter, reach_observations_file_name.
            It is expected to be in ASCII, comma-delimited format, and to contain two
            fields: gaged_reach_id (a string) and observed_value (a float).
        """
        self.gaged_reach_observations = {}
        f = open(reach_observations_file_name, 'r')
        # dictionary contain a list of grid-cells that constitute each gaged-reach
        for line in f.readlines():
            line_list = line.rstrip().split(',')
            gaged_reach_id, observed_value = tuple(line_list)
            observed_value = float(observed_value)
            self.gaged_reach_observations[gaged_reach_id] = observed_value
        f.close()

    def calc_gaged_reach_sim_flux(self, bc_reach_fluxes):
        """ Compute simulated flux for each gaged-reach. """
        self.gaged_reach_sim_fluxes = {}
        for gaged_reach_id in self.gaged_reach_list:
            for bc_type in self.gaged_reach_bc_reach_types[(gaged_reach_id)]:
                for bc_reach_id in self.gaged_reach_bc_reach_list[(gaged_reach_id, bc_type)]:
                    if self.gaged_reach_sim_fluxes.has_key((gaged_reach_id, bc_type)):
                        self.gaged_reach_sim_fluxes[(gaged_reach_id, bc_type)] += bc_reach_fluxes[(bc_type, bc_reach_id)]
                    else:
                        self.gaged_reach_sim_fluxes[(gaged_reach_id, bc_type)] = bc_reach_fluxes[(bc_type, bc_reach_id)]

    def output_gaged_reach_fluxes(self, sim_results_file_name):
        """ Output simulated and observed fluxes for each gaged reach. """
        
        bc_type_list = ['riv','drn','ghb']
        sim_results_file = self.init_obs_versus_sim_output_file(sim_results_file_name, bc_type_list)
        for gaged_reach_id in self.gaged_reach_list:
            observed_flux = self.gaged_reach_observations[gaged_reach_id]
            output_string = '{0}'.format(gaged_reach_id)
            total_sim_flux = 0
            for bc_type in bc_type_list:
                if self.gaged_reach_sim_fluxes.has_key((gaged_reach_id, bc_type)):
                    bc_flux = self.gaged_reach_sim_fluxes[(gaged_reach_id, bc_type)]
                else:
                    bc_flux = 0.
                total_sim_flux += bc_flux
                output_string += ',{0:.7e}'.format(bc_flux)
            residual = observed_flux - total_sim_flux
            if abs(observed_flux) > 1.e-7:
                residual_fraction = residual / observed_flux
            else:
                residual_fraction = -1.2345e25
            output_string += ',{0:.7e},{1:.7e},{2:.7e},{3:.7e}\n'.format(total_sim_flux, observed_flux, residual, residual_fraction)
            sim_results_file.write(output_string)
        sim_results_file.close()

    def init_obs_versus_sim_output_file(self, sim_results_file_name, bc_type_list):
        """ Open an output file with simulated versus observed flux data, and
            write a header with field names to that file.
        """
        sim_results_file = open(sim_results_file_name, 'w')
        output_string = 'gaged_reach_id'
        for bc_type in bc_type_list:
            output_string += ',{0}_sim_flux'.format(bc_type)
        output_string += ',sim_flux_total,obs_flux_total,residual,residual_fraction\n'
        sim_results_file.write(output_string)
        return(sim_results_file)
        
        
class ModflowListing:

    def __init__(self):
        #self.river_reach_list = []
        #self.river_reach_fluxes = {}
        #self.drain_list = []
        #self.drain_fluxes = {}
        #self.ghb_list = []
        #self.ghb_fluxes = {}
        bc_types = ['riv', 'drn', 'ghb']
        self.bc_reach_list = {}
        self.init_bc_reach_list(bc_types)
        self.bc_reach_fluxes = {}

    def init_bc_reach_list(self, bc_types):
        """ Initialize the dictionary of bc_reach_ids. """
        for bc_type in bc_types:
            self.bc_reach_list[bc_type] = []
        
    def openFiles(self):
        #inFileName = raw_input("Please enter name of MODFLOW input file for parsing: ")
        inFileName = 'NFSEG_SS_Y2001V1_141020updspr_2.lst'
        self.inFile = open(inFileName, 'r')
        outFileName = inFileName + '.out.csv'
        self.outFile = open(outFileName, 'w')

    def parseListingFile(self):
        while 1:
            line = self.inFile.readline()
            if not line:
                break
            self.checkForDrainFluxBlock(line)
            self.checkForRiverFluxBlock(line)
            self.checkForGHBFluxBlock(line)
        self.inFile.close()        

    def checkForGHBFluxBlock(self,line):
        if line.find("  HEAD DEP BOUNDS   PERIOD") == 0:
            stress_period = int(line.rstrip()[27:31])
            time_step = int(line.rstrip()[39:42])
            self.parseGHBFluxBlock() 

    def checkForDrainFluxBlock(self,line):
        if line.find("           DRAINS   PERIOD") == 0:
            stress_period = int(line.rstrip()[27:31])
            time_step = int(line.rstrip()[38:42])
            self.parseDrainFluxBlock()

    def checkForRiverFluxBlock(self,line):
        if line.find("    RIVER LEAKAGE   PERIOD") == 0:
            stress_period = int(line.rstrip()[27:31])
            time_step = int(line.rstrip()[38:42])
            self.parseRiverFluxBlock()

    def parseGHBFluxBlock(self):
        while 1:
            line = self.inFile.readline()
            if not line:
                raise ValueError, 'reached end of file while reading drain fluxes'
            elif line[:9] != ' BOUNDARY':
                break
            bc_reach_id = int(line.rstrip()[10:16])
            self.bc_reach_list['ghb'].append(bc_reach_id)
            flux = float(line.rstrip()[60:])
            self.bc_reach_fluxes[('ghb', bc_reach_id)] = flux

    def parseDrainFluxBlock(self):
        while 1:
            line = self.inFile.readline()
            if not line:
                raise ValueError, 'reached end of file while reading drain fluxes'
            elif line[:6] != ' DRAIN':
                break
            bc_reach_id = int(line.rstrip()[7:13])
            self.bc_reach_list['drn'].append(bc_reach_id)
            flux = float(line.rstrip()[60:])
            self.bc_reach_fluxes[('drn', bc_reach_id)] = flux

    def parseRiverFluxBlock(self):
        while 1:
            line = self.inFile.readline()
            if not line:
                raise ValueError, 'reached end of file while reading river fluxes'
            elif line[:6] != ' REACH':
                break
            bc_reach_id = int(line.rstrip()[7:13])
            self.bc_reach_list['riv'].append(bc_reach_id)
            flux = float(line.rstrip()[58:])
            self.bc_reach_fluxes[('riv', bc_reach_id)] = flux

    def close_files(self):
        #self.inFile.close()
        self.outFile.close()

def main():
    """ Compute simulated 'gaged-reach' fluxes by extract simulated drn, ghb,
        and riv fluxes from MODFLOW output listing. Compare them with observed
        values and output results to a file.
    """
    start_time = time.time()
    #shelf = shelve.open('sim_flux.shelf')

    modflow_input_format = 'free'
    output_file = open('gaged_reach_fluxes.csv', 'w')    

    # parse the MODFLOW listing file
    mf = ModflowListing()
    mf.openFiles()
    mf.parseListingFile()
    #shelf['bc_reach_fluxes'] = mf.bc_reach_fluxes
    
    # compute the simulated gaged-reach fluxes and compare with obs values
    gaged_reaches = GagedReaches('gaged_reach_definitions.csv')
    gaged_reaches.calc_gaged_reach_sim_flux(mf.bc_reach_fluxes)
    gaged_reaches.assign_observations_to_reaches('gaged_reach_observations.csv')
    gaged_reaches.output_gaged_reach_fluxes('gaged_reach_fluxes.csv')
    #shelf['gaged_reaches'] = gaged_reaches

    mf.close_files()
    #shelf.close()
    elapsed_time = time.time() - start_time
    print('\nAll Done!\n\nGaged-reach post-processing program executed in {0:0.1f} seconds.'.format(elapsed_time))

main()
