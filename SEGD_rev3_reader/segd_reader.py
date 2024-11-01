'''
`segd_reader.py`: Main SEG-D rev 3.0 reader sript.
Includes
- Class + functions to read and parse header/trace data from SEG-D rev 3.0 files.
- Function to return data as obspy stream.


MIT License

Copyright (c) 2024 Sacha Lapins

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
'''


import struct
import csv
import datetime
import numpy as np
from obspy import UTCDateTime, Trace, Stream
from pathlib import Path


''' 1. Main SEG-D rev 3.0 reader class '''

class SEG_D_Reader:
    ''' 
    Reads and parses various header/trace data and outputs as Python dictionary.

    Parameters:
    filepath: str, path to SEG-D rev 3.0 file to be read
    '''
    
    def __init__(self, filepath, segdreader_path = str(Path(__file__).parent), verbose=False):
        self.filepath = filepath
        self.segdreader_path = segdreader_path
        self.file = None
        self.verbose = verbose
    
    def open_file(self):
        """ Open SEG-D file in binary mode """
        try:
            self.file = open(self.filepath, 'rb')
        except IOError:
            print(f"Unable to open file {self.filepath}")
            return False
        return True

    def close_file(self):
        """ Close the opened SEG-D file """
        if self.file:
            self.file.close()

    def read_bcd(self, data, start_byte, end_byte, start_nibble, end_nibble):
        """
        Read Binary Coded Decimal (BCD) from the specified byte and nibble range.
        Uses 1-based indexing to align with fields given in SEG-D rev3.0 documentation.
        """
        # Extract the relevant bytes
        relevant_bytes = data[(start_byte-1):end_byte]

        # Convert bytes to a string of nibbles
        nibble_list = []
        for byte in relevant_bytes:
            # Convert each byte to its two nibbles (upper and lower nibbles)
            upper_nibble = (byte >> 4) & 0xF
            lower_nibble = byte & 0xF
            if upper_nibble < 10:
                nibble_list.append(str(upper_nibble))
            else:
                nibble_list.append(hex(upper_nibble)[2:])
            if lower_nibble < 10:
                nibble_list.append(str(lower_nibble))
            else:
                nibble_list.append(hex(lower_nibble)[2:])

        # Total nibbles in the relevant bytes
        total_nibbles = len(nibble_list)

        # Adjust start_nibble and end_nibble relative to the start of the extracted bytes
        nibble_offset = (start_byte - 1) * 2  # Since each byte gives two nibbles

        adjusted_start_nibble = start_nibble - nibble_offset - 1  # 1-based to 0-based index
        adjusted_end_nibble = end_nibble - nibble_offset

        # Extract the nibbles in the specified range
        nibble_string = ''.join(str(nibble_list[i]) for i in range(adjusted_start_nibble, adjusted_end_nibble))

        if nibble_string:
            # return int(nibble_string)  # Convert to integer if valid
            return nibble_string
        else:
            return 0  # Return 0 if the nibble string is empty

    def read_unsigned_binary(self, data, start_byte, end_byte):
        """ Read unsigned binary integer from the specified byte range (1-based indexing) """
        relevant_bytes = data[(start_byte-1):end_byte]
        return int.from_bytes(relevant_bytes, byteorder='big')

    def read_4bit_unsigned_binary(self, data, byte_pos, nibble):
        """
        Read a 4-bit unsigned binary number (nibble) from the specified byte position (1-based indexing).
        
        Parameters:
        - data: the binary data
        - byte_pos: the byte position (1-based index)
        - nibble: which nibble to read ('upper' or 'lower')
        
        Returns:
        - The 4-bit unsigned binary number.
        """
        byte_value = data[byte_pos - 1]  # Get the byte (adjust for 1-based index)
        
        if nibble == 'upper':
            return (byte_value >> 4) & 0x0F  # Upper nibble (first 4 bits)
        elif nibble == 'lower':
            return byte_value & 0x0F  # Lower nibble (last 4 bits)
        else:
            raise ValueError("Nibble must be 'upper' or 'lower'")

    def read_fractional_unsigned_binary(self, data, start_byte, end_byte):
        relevant_bytes = data[(start_byte-1):end_byte]
        fractional = int.from_bytes(relevant_bytes, byteorder='big') / 256.0  # 2^8 = 256, so divide by 256 to get the fraction
        
        return fractional

    def read_signed_binary(self, data, start_byte, end_byte, endian="big"):

        # Extract relevant bytes from the data (adjusting for 1-based indexing)
        relevant_bytes = data[(start_byte-1):end_byte]

        # Determine the size of the range
        byte_range_size = end_byte - start_byte + 1

        if byte_range_size == 3:
            # Three-byte two's complement signed binary number
            relevant_bytes = data[(start_byte - 1):end_byte]  # Extract 3 bytes
            
            # Sign-extend the 3-byte value to 4 bytes for two's complement handling
            if relevant_bytes[0] & 0x80:
                # Negative number, sign extend with 0xFF
                padded_bytes = b'\xFF' + relevant_bytes
            else:
                # Positive number, sign extend with 0x00
                padded_bytes = b'\x00' + relevant_bytes
            
            # Convert to a signed 4-byte integer
            signed_value = int.from_bytes(padded_bytes, byteorder='big', signed=True)

        else:
            # Convert the extracted bytes into a signed integer (big-endian)
            signed_value = int.from_bytes(relevant_bytes, byteorder='big', signed=True)

        return signed_value

    def read_ieee_float(self, data, start_byte, end_byte):
        """
        Read a 4-byte IEEE 754 floating-point number (32-bit float).
        
        Parameters:
        - data: The binary data containing the floating-point number.
        - start_byte: The starting byte (1-based index).
        - end_byte: The ending byte (should be start_byte + 3 for a 4-byte field).
        
        Returns:
        - The floating-point number as a float.
        """
        # Ensure we're reading exactly 4 bytes
        if end_byte - start_byte != 3:
            raise ValueError("The IEEE float must span exactly 4 bytes.")
        
        # Extract relevant bytes (adjusting for 1-based indexing)
        relevant_bytes = data[(start_byte - 1):end_byte]

        # Check if the bytes represent the NaN case (0x7FFFFFFF)
        if relevant_bytes == b'\x7F\xFF\xFF\xFF':
            return float('nan')  # Return NaN if the bytes match the NaN pattern
        
        # Use struct to unpack the 4 bytes as a 32-bit IEEE float ('>f' = big-endian float)
        float_value = struct.unpack('>f', relevant_bytes)[0]
        
        return float_value

    def read_ieee_double_float(self, data, start_byte, end_byte):
        """
        Read an 8-byte IEEE 754 floating-point number (64-bit float).
        
        Parameters:
        - data: The binary data containing the floating-point number.
        - start_byte: The starting byte (1-based index).
        - end_byte: The ending byte (should be start_byte + 7 for an 8-byte field).
        
        Returns:
        - The floating-point number as a float.
        """
        # Ensure we're reading exactly 8 bytes
        if end_byte - start_byte != 7:
            raise ValueError("The IEEE double float must span exactly 8 bytes.")
        
        # Extract relevant bytes (adjusting for 1-based indexing)
        relevant_bytes = data[(start_byte - 1):end_byte]
        
        # Use struct to unpack the 8 bytes as a 64-bit IEEE float ('>d' = big-endian double float)
        double_float_value = struct.unpack('>d', relevant_bytes)[0]
        
        return double_float_value
    
    def read_20byte_header(self, header_start):
        """ Generic function for reading 20 byte headers """
        if not self.file:
            return None

        self.file.seek(header_start) # Start at header_start
        header = self.file.read(20) # Read 20 bytes
        
        return header
    
    def read_32byte_header(self, header_start):
        """ Generic function for reading 32 byte headers """
        if not self.file:
            return None

        self.file.seek(header_start) # Start at header_start
        header = self.file.read(32) # Read 32 bytes
        
        return header

    def read_96byte_header(self, header_start):
        """ Generic function for reading 96 byte headers """
        if not self.file:
            return None

        self.file.seek(header_start) # Start at header_start
        header = self.file.read(96) # Read 96 bytes
        
        return header

    def read_header_block_type(self, header_start):
        """ Reads byte 32 from header to extract header block type code """
        if not self.file:
            return None

        self.file_seek(header_start + 31)
        header_block_byte = self.file.read(1)

        return self.read_unsigned_binary(header_block_byte, 1, 1)

    def read_segd_timestamp(self, data, start_byte, end_byte, start_nibble=None, end_nibble=None):
        """
        Read an 8-byte SEG-D timestamp (signed, big-endian integer) from the specified byte and nibble range.
        
        Parameters:
        - data: The binary data containing the SEG-D timestamp.
        - start_byte: The starting byte (1-based index).
        - end_byte: The ending byte (1-based index).
        - start_nibble: Optional, 'upper' or 'lower' to specify which nibble to start with in the first byte.
        - end_nibble: Optional, 'upper' or 'lower' to specify which nibble to end with in the last byte.
        
        Returns:
        - The SEG-D timestamp in microseconds as a signed integer.
        """
        # Extract relevant bytes from the data (adjusting for 1-based indexing)
        relevant_bytes = data[(start_byte-1):end_byte]
        
        # Handle starting nibble if specified (upper nibble = most significant, lower nibble = least significant)
        if start_nibble == 'lower':
            relevant_bytes[0] = relevant_bytes[0] & 0x0F  # Mask to keep only the lower nibble
        elif start_nibble == 'upper':
            relevant_bytes[0] = relevant_bytes[0] & 0xF0  # Mask to keep only the upper nibble
        
        # Handle ending nibble if specified
        if end_nibble == 'lower':
            relevant_bytes[-1] = relevant_bytes[-1] & 0x0F  # Mask to keep only the lower nibble
        elif end_nibble == 'upper':
            relevant_bytes[-1] = relevant_bytes[-1] & 0xF0  # Mask to keep only the upper nibble
        
        # Convert the extracted bytes into a signed 8-byte integer (big-endian)
        timestamp = int.from_bytes(relevant_bytes, byteorder='big', signed=True)
        
        return timestamp
    
    def calculate_leap_seconds(self, days_since_epoch):
        """Calculate leap seconds based on days since GPS epoch."""
        for days, leap_sec in reversed(LEAP_SECONDS):
            if days_since_epoch >= days:
                return leap_sec
        return 0
    
    def segd_timestamp_to_UCTDateTime(self, gps_timestamp):
        """Convert SEG-D GPS timestamp to UTC."""
        gps_epoch = datetime.datetime(1980, 1, 6, 0, 0, 0)
        
        # Calculate the total seconds from microseconds
        total_seconds = gps_timestamp / 1e6
        
        # Calculate days and remaining seconds since GPS epoch
        days_since_epoch = int(total_seconds // 86400)
        remaining_seconds = total_seconds % 86400
        
        # Get leap seconds correction
        leap_seconds = self.calculate_leap_seconds(days_since_epoch)
        
        # Correct the timestamp by subtracting leap seconds
        utc_seconds = total_seconds - leap_seconds
        
        # Convert back to datetime format
        utc_time = gps_epoch + datetime.timedelta(seconds=utc_seconds)
        
        return UTCDateTime(utc_time)

    def read_ascii(self, data, start_byte, end_byte):
        """
        Reads ASCII text, left-justified and space-padded.
        
        Parameters:
        - data: The binary data containing the serial number.
        - start_byte: The starting byte (1-based index).
        - end_byte: The ending byte.
        
        Returns:
        - String containing text.
        """    
        # Extract relevant bytes (adjusting for 1-based indexing)
        relevant_bytes = data[(start_byte - 1):end_byte]
        
        # Decode the bytes as ASCII and strip trailing spaces (padded with 0x20)
        out_string = relevant_bytes.decode('ascii').rstrip(' ')
        
        return out_string
    
    def read_header_spec_from_csv(self, spec_filepath):
        """ Reads header specification from CSV file and returns it as a list of dictionaries """
        header_spec = []
        with open(spec_filepath, 'r') as csvfile:
            csv_reader = csv.DictReader(csvfile)
            for row in csv_reader:
                header_spec.append({
                    'name': row['name'],
                    'start_byte': int(row['startByte']),
                    'end_byte': int(row['endByte']),
                    'start_nibble': int(row['startNibble']),
                    'end_nibble': int(row['endNibble']),
                    'format': row['format']
                })
        return header_spec
    
    def read_header_value_from_csv(self, value_filepath, value):  
        """ Reads a CSV file containing code-description pairs and returns the description for the provided code. """        
        with open(value_filepath, mode='r') as csvfile:
            csv_reader = csv.DictReader(csvfile)

            for row in csv_reader:
                if row['code'] == value:
                    return value + ": " + str(row['description'])  # Return the matched description
                
        return value + ": Unknown code"  # If the code is not found, return 'Unknown code'
    
    def parse_header_fields_from_spec(self, header, header_spec):
        """ Use dictionaries from read_header_spec_from_csv to parse header fields """
        # Initialize a dictionary to store parsed header fields
        header_fields = {}

        # Iterate over the header specification and parse fields
        for spec in header_spec:
            field_name = spec['name']
            start_byte = spec['start_byte']
            end_byte = spec['end_byte']
            start_nibble = spec['start_nibble']
            end_nibble = spec['end_nibble']
            data_format = spec['format']

            # Parse based on the format (bcd or ubin)
            if data_format == 'bcd':
                value = self.read_bcd(header, start_byte, end_byte, start_nibble, end_nibble)
            elif data_format == 'ubin':
                if start_nibble == end_nibble:
                    if start_nibble == (2 * start_byte) - 1:
                        value = self.read_4bit_unsigned_binary(header, start_byte, 'upper')
                    else:
                        value = self.read_4bit_unsigned_binary(header, start_byte, 'lower')
                    if field_name == "generalHeaderBlocks":
                        # if value > 14:
                            # value = "FF"
                        value = hex(value)[2:]
                else:
                    value = self.read_unsigned_binary(header, start_byte, end_byte)
                    if field_name == "baseScanInterval":
                        value = value * 1/16
                        # if value == 255/16:
                        #     value = "ff"
                    # if field_name == "dominantSamplingInterval":
                        # value = value * 1e-6
                # if field_name == "headerBlockType":
                #     value = hex(value)
            elif data_format == 'sbin':
                value = self.read_signed_binary(header, start_byte, end_byte)
            elif data_format == 'fraction':
                value = self.read_fractional_unsigned_binary(header, start_byte, end_byte)
            elif data_format == "timestamp":
                value = self.read_segd_timestamp(header, start_byte, end_byte)
            elif data_format == "ieee":
                value = self.read_ieee_float(header,start_byte,end_byte)
            elif data_format == "double":
                value = self.read_ieee_double_float(header,start_byte,end_byte)
            elif (data_format == "serial") or (data_format == "ascii"):
                value = self.read_ascii(header,start_byte,end_byte)
            elif data_format == "hex":
                value = hex(self.read_unsigned_binary(header,start_byte,end_byte))
            else:
                raise ValueError(f"Unsupported format {data_format}")

            # Store the parsed value in the header_fields dictionary
            header_fields[field_name] = value
        
        return header_fields

    def read_storage_unit_label(self):
        """ Reads the first 128 bytes of ASCII storage unit label """
        if not self.file:
            return None

        # Read the first 128 bytes which is the storage unit label
        self.file.seek(0)
        storage_unit_label = self.file.read(128).decode('ascii', errors='ignore').strip()
        
        return storage_unit_label        
    
    def read_general_header1(self, after_storage_unit_label=False):
        """ Reads General Header #1 """
        if not self.file:
            return None

        # Read General Header #1
        if after_storage_unit_label:
            general_header_start = 128  # Start after the 128-byte storage unit label
        else:
            general_header_start = 0
        general_header = self.read_32byte_header(general_header_start)

        # Parse fields according to CSV specification
        general_header_spec = self.read_header_spec_from_csv(self.segdreader_path + "/segd_csv_headers/generalHeader1.csv")
        header_fields = self.parse_header_fields_from_spec(general_header, general_header_spec)

        if self.verbose:
            # Print parsed fields for verification
            for field, value in header_fields.items():
                print(f"{field}: {value}")
        
        return header_fields
    
    def read_general_header2(self, include_storage_unit_label=False):
        """ Reads General Header #2 """
        if not self.file:
            return None

        # Read General Header #2
        if include_storage_unit_label:
            general_header_start = 128 + 32  # Start after the 128-byte storage unit label
        else:
            general_header_start = 32
        general_header = self.read_32byte_header(general_header_start)

        # Parse fields according to CSV specification
        general_header_spec = self.read_header_spec_from_csv(self.segdreader_path + "/segd_csv_headers/generalHeader2.csv")
        header_fields = self.parse_header_fields_from_spec(general_header, general_header_spec)

        if self.verbose:
            # Print parsed fields for verification
            for field, value in header_fields.items():
                print(f"{field}: {value}")
        
        return header_fields

    def read_general_header3(self, include_storage_unit_label=False):
        """ Reads General Header #3 """
        if not self.file:
            return None

        # Read General Header #3
        if include_storage_unit_label:
            general_header_start = 128 + 64  # Start after the 128-byte storage unit label
        else:
            general_header_start = 64
        general_header = self.read_32byte_header(general_header_start)

        # Parse fields according to CSV specification
        general_header_spec = self.read_header_spec_from_csv(self.segdreader_path + "/segd_csv_headers/generalHeader3.csv")
        header_fields = self.parse_header_fields_from_spec(general_header, general_header_spec)

        if self.verbose:
            # Print parsed fields for verification
            for field, value in header_fields.items():
                print(f"{field}: {value}")
        
        return header_fields

    def read_scan_type_header(self, start_byte):
        """ Reads Scan Type Header """
        if not self.file:
            return None

        scan_header = self.read_96byte_header(start_byte)

        # Parse fields according to CSV specification
        scan_header_spec = self.read_header_spec_from_csv(self.segdreader_path + "/segd_csv_headers/scanTypeHeader.csv")
        header_fields = self.parse_header_fields_from_spec(scan_header, scan_header_spec)

        if self.verbose:
            # Print parsed fields for verification
            for field, value in header_fields.items():
                print(f"{field}: {value}")
        
        return header_fields

    def read_demux_trace_header(self, start_byte):
        """ Reads Demux trace header """
        if not self.file:
            return None

        demux_header = self.read_20byte_header(start_byte)

        # Parse fields according to CSV specification
        demux_header_spec = self.read_header_spec_from_csv(self.segdreader_path + "/segd_csv_headers/demuxTraceHeader.csv")
        header_fields = self.parse_header_fields_from_spec(demux_header, demux_header_spec)

        if self.verbose:
            # Print parsed fields for verification
            for field, value in header_fields.items():
                print(f"{field}: {value}")
        
        return header_fields

    def read_trace_header_extension1(self, start_byte):
        """ Reads Trace Header Extension #1 """
        if not self.file:
            return None

        trace_header = self.read_32byte_header(start_byte)

        # Parse fields according to CSV specification
        trace_header_spec = self.read_header_spec_from_csv(self.segdreader_path + "/segd_csv_headers/traceHeaderExtension.csv")
        header_fields = self.parse_header_fields_from_spec(trace_header, trace_header_spec)

        if self.verbose:
            # Print parsed fields for verification
            for field, value in header_fields.items():
                print(f"{field}: {value}")
        
        return header_fields

    def read_other_header(self, start_byte):
        """ Reads other / optional header extensions given in SEG-D rev 3.0 documentation """
        if not self.file:
            return None

        trace_header = self.read_32byte_header(start_byte)

        # Read byte32 and check header block type:
        byte32 = self.read_unsigned_binary(trace_header, 32, 32)

        # Use header code to identify header type and csv file
        if byte32 == int('0x10',16):
            trace_header_csv_file = 'generalHeader4.csv'
        elif byte32 == int('0x11',16):
            trace_header_csv_file = 'generalHeader5.csv'
        elif byte32 == int('0x12',16):
            trace_header_csv_file = 'generalHeader6.csv'
        elif byte32 == int('0x30',16):
            trace_header = self.read_96byte_header(start_byte) # hex 30 is code for 96 byte scan type header (channel set description)
            trace_header_csv_file = 'scanTypeHeader.csv'
        elif byte32 == int('0x31',16):
            # This is part of above channel set description, skip
            return None
        elif byte32 == int('0x32',16):
            # This is part of above channel set description, skip
            return None
        elif byte32 == int('0x41',16):
            trace_header_csv_file = 'sensorInfoHeader.csv'
        elif byte32 == int('0x42',16):
            trace_header_csv_file = 'timestampHeader.csv'
        elif byte32 == int('0x43',16):
            trace_header_csv_file = 'sensorCalibrationHeader.csv'
        elif byte32 == int('0x44',16):
            trace_header_csv_file = 'timeDriftHeader.csv'
        elif byte32 == int('0x50',16):
            trace_header = self.read_96byte_header(start_byte) # hex 50 is code for 96 byte position blocks header
            trace_header_csv_file = 'positionBlocks.csv'
        elif byte32 == int('0x51',16):
            # This is part of above position blocks header, skip
            return None
        elif byte32 == int('0x52',16):
            # This is part of above position blocks header, skip
            return None
        elif byte32 == int('0x55',16):
            trace_header_csv_file = 'coordRefSystem.csv'
        elif byte32 == int('0x56',16):
            trace_header_csv_file = 'relativePositionBlock.csv'
        elif byte32 == int('0x60',16):
            trace_header_csv_file = 'orientationHeader.csv'
        elif byte32 == int('0x61',16):
            trace_header_csv_file = 'measurementBlock.csv'
        elif byte32 == int('0xd6',16):
            trace_header_csv_file = 'userD6.csv'
        else:
            if self.verbose:
                print("Header type unrecognised, skipping. Header block type: " + hex(byte32)[2:])
            return None
            
        # Parse fields according to CSV specification
        trace_header_spec = self.read_header_spec_from_csv(self.segdreader_path + "/segd_csv_headers/" + trace_header_csv_file)
        header_fields = self.parse_header_fields_from_spec(trace_header, trace_header_spec)

        if self.verbose:
            # Print parsed fields for verification
            for field, value in header_fields.items():
                print(f"{field}: {value}")
        
        return header_fields


    def read_32bit_IEEE_trace_data(self, start_byte, num_samples):
        """
        Read trace data from a SEG-D file.
        
        Parameters:
        - data: Binary data containing the trace samples.
        - start_byte: The starting byte of the trace data (1-based index).
        - num_samples: Number of samples in the trace.
        
        Returns:
        - A list of decoded trace samples.
        """        

        if not self.file:
            return None

        self.file.seek(start_byte) # Start of trace data
        total_trace_bytes = num_samples * 4 # Total size of trace data in bytes
        relevant_bytes = self.file.read(total_trace_bytes) # Read required no. of bytes for trace
        
        # Unpack the bytes into a list of samples
        fmt = f'>{num_samples}f'  # Big-endian 32-bit IEEE floats
        trace_samples = struct.unpack(fmt, relevant_bytes)
        
        return trace_samples

    def get_data_plus_key_headers(self):
        """
        Stores key header from array in dictionary
        """
        if not self.file:
            try:
                self.open_file()
            except:
                print("Can't open file.")

        out_dict = {}
        out_dict['file_header'] = {}

        if self.verbose:
            print("*** Reading file headers ***")
        # Read first three general headers and extract key info
        gen_head1 = self.read_general_header1()
        gen_head2 = self.read_general_header2()
        gen_head3 = self.read_general_header3()

        fileNumber = gen_head1['fileNumber']
        if fileNumber == "ffff":
            fileNumber = gen_head2['expandedFileNumber']

        duration = gen_head1['recordLength']
        if 'f' in duration:
            duration = gen_head2['extendedRecordLength']
        
        samplingInterval = gen_head1['baseScanInterval']
        if samplingInterval == 255/16:
            samplingInterval = gen_head2['dominantSamplingInterval']
        
        out_dict['file_header']['fileNumber'] = fileNumber
        out_dict['file_header']['SEGDVersion'] = str(gen_head2['majorSEGDrevisionNumber']) + "." + str(gen_head2['minorSEGDrevisionNumber'])
        out_dict['file_header']['year'] = gen_head1['year']
        out_dict['file_header']['jday'] = gen_head1['day']
        out_dict['file_header']['hour'] = gen_head1['hour']
        out_dict['file_header']['minute'] = gen_head1['minute']
        out_dict['file_header']['second'] = gen_head1['second']
        out_dict['file_header']['duration_microsec'] = duration
        out_dict['file_header']['dominant_sampling_interval_microsec'] = samplingInterval
        out_dict['file_header']['record_timezero_utc'] = self.segd_timestamp_to_UCTDateTime(gen_head3['timeZero'])

        headerSize = gen_head3['headerSize']
        out_dict['file_header']['headerSize_bytes'] = headerSize

        # Loop through remaining general header blocks while i < headerSize
        i = 3*32 # Start after first three general headers:
        num_channel_sets = 0
        while i < headerSize:
            tmp = self.read_other_header(i) # Returns None if unrecognised header, otherwise dict
            if type(tmp) is dict:
                if tmp['headerBlockType'] == 48:
                    out_dict['channelSet_' + str(num_channel_sets + 1)] = {}
                    out_dict['channelSet_' + str(num_channel_sets + 1)]['description'] = tmp
                    out_dict['channelSet_' + str(num_channel_sets + 1)]['traceData'] = {}
                    num_channel_sets += 1
            i += 32

        # Then loop through channel sets (trace data)
        if self.verbose:
            print("*** Reading trace headers and data ***")
        start_byte = headerSize # start place in file
        for channel_set in range(num_channel_sets):
            num_trace_headers = out_dict['channelSet_' + str(channel_set + 1)]['description']['numberTraceHeaderExtensions']
            num_traces = out_dict['channelSet_' + str(channel_set + 1)]['description']['numberOfChannelsThisSet']
            num_samples_per_trace = out_dict['channelSet_' + str(channel_set + 1)]['description']['samplesPerChannel']
            # Loop through traces
            for trace in range(num_traces):
                if self.verbose:
                    print("*** Reading trace number " + str(trace + 1) + " of " + str(num_traces) + " (channel set #" + str(channel_set + 1) + " of " + str(num_channel_sets) + ") ***")
                tmp = self.read_demux_trace_header(start_byte)
                start_byte += 20
                # Loop through trace headers:
                for trace_header in range(num_trace_headers):
                    if self.verbose:
                        print("*** Reading trace header " + str(trace_header + 1) + " of " + str(num_trace_headers) + " ***")
                    if trace_header == 0:
                        tmp = self.read_trace_header_extension1(start_byte)
                        line_num = tmp['receiverLineNumber']
                        if line_num < 0:
                            line_num = tmp['extendedReceiverLineNumberInteger']
                        line_num = 'line_' + str(line_num)
                        point_num = tmp['receiverPointNumber']
                        if point_num < 0:
                            point_num = tmp['extendedReceiverPointNumberInteger']
                        point_num = 'point_' + str(point_num)
                        if line_num not in out_dict['channelSet_' + str(channel_set + 1)]['traceData']:
                            out_dict['channelSet_' + str(channel_set + 1)]['traceData'][line_num] = {}
                        if point_num not in out_dict['channelSet_' + str(channel_set + 1)]['traceData'][line_num]:
                            out_dict['channelSet_' + str(channel_set + 1)]['traceData'][line_num][point_num] = {}
                        out_dict['channelSet_' + str(channel_set + 1)]['traceData'][line_num][point_num]['trace_header'] = {}
                        out_dict['channelSet_' + str(channel_set + 1)]['traceData'][line_num][point_num]['trace_header']['sensorType_code'] = tmp['sensorType']
                        out_dict['channelSet_' + str(channel_set + 1)]['traceData'][line_num][point_num]['trace_header']['physicalUnit_code'] = tmp['physicalUnit']
                        out_dict['channelSet_' + str(channel_set + 1)]['traceData'][line_num][point_num]['trace_header']['numSamplesPerTrace'] = tmp['numSamplesPerTrace']
                    else:
                        tmp = self.read_other_header(start_byte)
                        if type(tmp) is dict:
                            if tmp['headerBlockType'] == 66:
                                out_dict['channelSet_' + str(channel_set + 1)]['traceData'][line_num][point_num]['trace_header']['timeZero_utc'] = self.segd_timestamp_to_UCTDateTime(tmp['timeZero'])
                            elif tmp['headerBlockType'] == 65:
                                out_dict['channelSet_' + str(channel_set + 1)]['traceData'][line_num][point_num]['trace_header']['sensorSensitivity'] = tmp['sensorSensitivity']
                                out_dict['channelSet_' + str(channel_set + 1)]['traceData'][line_num][point_num]['trace_header']['serialNumber'] = tmp['serialNumber']
                            elif tmp['headerBlockType'] == 214:
                                out_dict['channelSet_' + str(channel_set + 1)]['traceData'][line_num][point_num]['trace_header']['latitude'] = tmp['latitude']
                                out_dict['channelSet_' + str(channel_set + 1)]['traceData'][line_num][point_num]['trace_header']['longitude'] = tmp['longitude']
                                out_dict['channelSet_' + str(channel_set + 1)]['traceData'][line_num][point_num]['trace_header']['elevation'] = tmp['elevation']
                    start_byte += 32
                # Get trace data:
                out_dict['channelSet_' + str(channel_set + 1)]['traceData'][line_num][point_num]['trace_data'] = self.read_32bit_IEEE_trace_data(start_byte, out_dict['channelSet_' + str(channel_set + 1)]['traceData'][line_num][point_num]['trace_header']['numSamplesPerTrace'])
                start_byte += out_dict['channelSet_' + str(channel_set + 1)]['traceData'][line_num][point_num]['trace_header']['numSamplesPerTrace'] * 4
                
        return out_dict


def SEG_D_to_stream(filelist, convert_to_int = True, serial_to_station_name_dict = None, network_code = 'AA', remove_gaps = False):
    '''
    Reads a Sercel SEG-D revision 3.0 file and returns an obspy stream
    '''

    # Initiate obspy stream for all data in filelist
    st = Stream()
    
    if type(filelist) is str:
        filelist = [filelist]

    for filepath in filelist:
        # Initiate SEG_D_Reader:
        reader = SEG_D_Reader(filepath)
    
        # Open file:
        reader.open_file()
    
        # Get key header info and trace data:
        data = reader.get_data_plus_key_headers()
    
        channel_sets = [key for key in list(data.keys()) if 'channelSet' in key]
        channel_sets.sort()
        line_nums = [list(data[chan_set]['traceData'].keys()) for chan_set in channel_sets]
    
        for chan_set, line_num_list in zip(channel_sets, line_nums):
            for line_num in line_num_list:
                for trace_num in list(data[chan_set]['traceData'][line_num].keys()):
                    trace_h, trace_d = data[chan_set]['traceData'][line_num][trace_num]['trace_header'], np.array(data[chan_set]['traceData'][line_num][trace_num]['trace_data'], dtype=np.float32)

                    if remove_gaps:
                        if trace_h['timeZero_utc'] == UTCDateTime("1980-01-06T00:00:00.000000Z"):
                            continue
                
                    # check if all traces can be converted to int
                    convert_to_int = convert_to_int and np.all(np.mod(trace_d, 1) == 0)
    
                    tr = Trace(trace_d)
                    if convert_to_int:
                        tr.data = tr.data.astype(np.int32)
                    if serial_to_station_name_dict is not None:
                        if 'network' in serial_to_station_name_dict[trace_h['serialNumber']]:
                            tr.stats.network = serial_to_station_name_dict[trace_h['serialNumber']]['network_code']
                        else:
                            tr.stats.network = network_code
                        tr.stats.station = serial_to_station_name_dict[trace_h['serialNumber']]['station_code']
                    else:
                        tr.stats.network = network_code
                        tr.stats.station = str(trace_h['serialNumber'])
    
                    sample_rate = data[chan_set]['description']['samplingInterval'] / 1e6 # Given in microseconds
                    tr.stats.delta = sample_rate
                    
                    if (1./sample_rate) >= 1000:
                        tr.stats.channel = 'G'
                    elif (1./sample_rate) >= 250:
                        tr.stats.channel = 'D'
                    elif (1./sample_rate) >= 80:
                        tr.stats.channel = 'E'
                    elif (1./sample_rate) >= 10:
                        tr.stats.channel = 'S'
    
                    instrument_code = trace_h['sensorType_code']
                    if serial_to_station_name_dict is not None:
                        if len(SENSOR_CODE[instrument_code]) > 1:
                            tr.stats.channel += SENSOR_CODE[instrument_code][:-1]
                        tr.stats.channel += serial_to_station_name_dict[trace_h['serialNumber']]['component']
                    else:
                        tr.stats.channel += SENSOR_CODE[instrument_code]
                
                    # tr.stats.starttime = trace_h['timeZero_utc']
                    tr.stats.starttime = data['file_header']['record_timezero_utc']
                    tr.stats.segd = {}
                    tr.stats.segd['serialNumber'] = str(trace_h['serialNumber'])
                    tr.stats.segd.update(data['file_header'])
                    tr.stats.segd['gainControl'] = CHANNEL_GAIN_CONTROL_CODE[data[chan_set]['description']['channelGainControl']]
                    tr.stats.segd['aliasFilterFrequency'] = data[chan_set]['description']['aliasFilterFrequency']
                    tr.stats.segd['lowCutFilterFrequency'] = data[chan_set]['description']['lowCutFilterFrequency']
                    tr.stats.segd['aliasFilterSlope'] = data[chan_set]['description']['aliasFilterSlope']
                    tr.stats.segd['lowCutFilterSlope'] = data[chan_set]['description']['lowCutFilterSlope']
                    tr.stats.segd['notchFrequency'] = data[chan_set]['description']['notchFrequency']
                    tr.stats.segd['secondNotchFrequency'] = data[chan_set]['description']['secondNotchFrequency']
                    tr.stats.segd['thirdNotchFrequency'] = data[chan_set]['description']['thirdNotchFrequency']
                    tr.stats.segd['filterPhase'] = FILTER_PHASE_CODE[data[chan_set]['description']['filterPhase']]
                    tr.stats.segd['filterDelaySecs'] = data[chan_set]['description']['filterDelay'] / 1e6
                    tr.stats.segd['DSM'] = data[chan_set]['description']['descaleMultiplier']
                    tr.stats.segd['sensitivity'] = trace_h['sensorSensitivity']
                    tr.stats.segd['physicalUnit'] = PHYSICAL_UNIT_CODE[trace_h['physicalUnit_code']]
                    tr.stats.segd['latitude'] = trace_h['latitude']
                    tr.stats.segd['longitude'] = trace_h['longitude']
                    tr.stats.segd['elevation'] = trace_h['elevation']
                    st.append(tr)
                    
        reader.close_file()
       
    st.merge() # Merge consecutive files

    return st



SENSOR_CODE = {
    0: '',    # not defined
    1: 'DH',  # hydrophone
    2: 'HZ',  # geophone, vertical
    3: 'H1',  # geophone, horizontal, in-line
    4: 'H2',  # geophone, horizontal, crossline
    5: 'H3',  # geophone, horizontal, other
    6: 'NZ',  # accelerometer, vertical
    7: 'N1',  # accelerometer, horizontal, in-line
    8: 'N2',  # accelerometer, horizontal, crossline
    9: 'N3'   # accelerometer, horizontal, other
}

PHYSICAL_UNIT_CODE = {
    0: '',         # not defined
    1: 'mbar',     # millibar
    2: 'bar',      # bar
    3: 'mm/s',     # millimeter/second
    4: 'm/s',      # meter/second
    5: 'mm/s/s',   # millimeter/second/second
    6: 'm/s/s',    # meter/second/second
    7: 'N',        # Newton
    8: 'K',        # Kelvin
    9: 'Hz',       # Hertz
    10: 's',       # second
    11: 'T',       # Tesla
    12: 'V/m',     # Volt/meter
    13: 'Vm',      # Volt meter
    14: 'A/m',     # Ampere/meter
    15: 'V',       # Volt
    16: 'A',       # Ampere
    17: 'Rad',     # Radians 
    18: 'Pa',      # Pascal
    19: 'mu_Pa',   # Micropascal
    20: 'mm',      # millimeter
    21: 'm'        # meter
}

FILTER_PHASE_CODE = {
    0: 'Unknown',
    1: 'Minimum',
    2: 'Linear',
    3: 'Zero',
    4: 'Mixed',
    5: 'Maximum'
}

CHANNEL_GAIN_CONTROL_CODE = {
    1: 'Individual AGC',
    2: 'Ganged AGC',
    3: 'Fixed gain',
    4: 'Programmed gain',
    8: 'Binary gain control',
    9: 'IFP gain control'
}

# Leap seconds data (days since GPS epoch, leap seconds added)
LEAP_SECONDS = [
    (542, 1),    # 1 July 1981
    (907, 2),    # 1 July 1982
    (1272, 3),   # 1 July 1983
    (2003, 4),   # 1 July 1985
    (2917, 5),   # 1 January 1988
    (3648, 6),   # 1 January 1990
    (4013, 7),   # 1 January 1991
    (4560, 8),   # 1 July 1992
    (4925, 9),   # 1 July 1993
    (5290, 10),  # 1 July 1994
    (5839, 11),  # 1 January 1996
    (6386, 12),  # 1 July 1997
    (6935, 13),  # 1 January 1999
    (9492, 14),  # 1 January 2006
    (10588, 15), # 1 January 2009
    (11865, 16), # 1 July 2012
    (12984, 17), # 1 July 2015
    (13244, 18)  # 1 January 2017
]