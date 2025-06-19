'''
`segd_reader.py`: Main SEG-D reader sript.

Currently supports the following SEG-D versions:
- rev 3.0
- rev 2.1

Tested on data from the following manufacturers:
- Sercel: rev 3.0
- Stryde: rev 3.0
- SmartSolo: rev 2.1

Includes
- Class + functions to read and parse header/trace data from SEG-D files.
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
import json
import datetime
import numpy as np
from obspy import UTCDateTime, Trace, Stream
from pathlib import Path
import importlib.resources as pkg_resources


''' 1. Main SEG-D reader class '''

class SEG_D_Reader:
    ''' 
    Reads and parses various header/trace data and outputs as Python dictionary.

    Parameters:
    filepath: str, path to SEG-D file to be read
    '''
    
    def __init__(self, filepath, verbose=False):
        self.filepath = filepath
        self.file = None
        self.major_version = None
        self.minor_version = None
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

    def load_header_specs(self, json_path):
        with open(json_path, 'r') as f:
            self.header_specs = json.load(f)

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
    
    def segd_timestamp_to_UTCDateTime(self, gps_timestamp):
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

    def read_header_spec_from_json(self, header_specs, key):
        """
        Retrieves the header spec dictionary for a given header key from header JSON file.
        
        Parameters:
            header_specs: dict loaded from header_specs.json
            key: string like '0x10' or '_generalHeader1'
    
        Returns:
            List of field spec dictionaries
        """
        try:
            return header_specs[key]["fields"]
        except KeyError:
            raise ValueError(f"Header spec for key '{key}' not found in JSON.")    
    

    def parse_header_fields_from_spec(self, header, header_spec):
        """
        Parses the SEG-D header based on JSON header field specs.
        
        Parameters:
            header: raw bytes of the header
            header_spec: list of field dictionaries from JSON
    
        Returns:
            dict of header fields
        """
        header_fields = {}
    
        for spec in header_spec:
            field_name = spec['name']
            start_byte = spec['startByte']
            end_byte = spec['endByte']
            start_nibble = spec['startNibble']
            end_nibble = spec['endNibble']
            data_format = spec['format']
    
            if data_format == 'bcd':
                value = self.read_bcd(header, start_byte, end_byte, start_nibble, end_nibble)
    
            elif data_format == 'ubin':
                if start_nibble == end_nibble:
                    if start_nibble == (2 * start_byte) - 1:
                        value = self.read_4bit_unsigned_binary(header, start_byte, 'upper')
                    else:
                        value = self.read_4bit_unsigned_binary(header, start_byte, 'lower')
                    if field_name == "generalHeaderBlocks":
                        value = hex(value)[2:]
                else:
                    value = self.read_unsigned_binary(header, start_byte, end_byte)
                    if field_name == "baseScanInterval":
                        value = value * (1 / 16)
    
            elif data_format == 'sbin':
                value = self.read_signed_binary(header, start_byte, end_byte)
    
            elif data_format == 'fraction':
                value = self.read_fractional_unsigned_binary(header, start_byte, end_byte)
    
            elif data_format == "timestamp":
                value = self.read_segd_timestamp(header, start_byte, end_byte)
    
            elif data_format == "ieee":
                value = self.read_ieee_float(header, start_byte, end_byte)
    
            elif data_format == "double":
                value = self.read_ieee_double_float(header, start_byte, end_byte)
    
            elif data_format in ("serial", "ascii"):
                value = self.read_ascii(header, start_byte, end_byte)
    
            elif data_format == "hex":
                value = hex(self.read_unsigned_binary(header, start_byte, end_byte))
    
            else:
                raise ValueError(f"Unsupported format '{data_format}' in field '{field_name}'")
    
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

    
    def read_and_parse_header(self, start_byte, byte_size, json_key):
        """Generalized SEG-D header reader using JSON spec."""
        if not self.file:
            return None
    
        read_method = {
            20: self.read_20byte_header,
            32: self.read_32byte_header,
            96: self.read_96byte_header
        }.get(byte_size)
    
        if not read_method:
            raise ValueError(f"Unsupported header size: {byte_size} bytes")
    
        header = read_method(start_byte)
    
        header_spec = self.read_header_spec_from_json(self.header_specs, json_key)
        header_fields = self.parse_header_fields_from_spec(header, header_spec)
    
        if self.verbose:
            print(f"Header type: {json_key}")
            for field, value in header_fields.items():
                print(f"{field}: {value}")
    
        return header_fields

    
    def read_general_header1(self, include_storage_unit_label=False):
        """ Reads General Header #1 """
        offset = 128 if include_storage_unit_label else 0        
        return self.read_and_parse_header(offset, 32, "_generalHeader1")


    def read_general_header2(self, include_storage_unit_label=False):
        """ Reads General Header #2 """
        offset = 128 + 32 if include_storage_unit_label else 32
        return self.read_and_parse_header(offset, 32, "_generalHeader2")


    def read_general_header3(self, include_storage_unit_label=False):
        """ Reads General Header #3 """
        if int(self.major_version) < 3:
            gen_head1 = self.read_general_header1()
            if int(gen_head1.get("generalHeaderBlocks", 0), 16) < 2:
                return None
            return self.read_and_parse_header(64 + (128 if include_storage_unit_label else 0), 32, "_generalHeaderN")  # Update "0xNN" if you define one
        return self.read_and_parse_header(64 + (128 if include_storage_unit_label else 0), 32, "0x03")



    def read_general_headerN(self, start_byte):
        """ Reads General Header N (only applies to rev 2.1, and possibly below - not checked earlier versions) """
        if int(self.major_version) > 2:
            return None
        return self.read_and_parse_header(start_byte, 32, "_generalHeaderN")  # Assign a key for GHN in the v2.1 JSON


    def read_scan_type_header(self, start_byte):
        """ Reads Scan Type Header """
        if int(self.major_version) > 2:
            return self.read_and_parse_header(start_byte, 96, "0x30")
        else:
            return self.read_and_parse_header(start_byte, 32, "0x30")

    
    def read_demux_trace_header(self, start_byte):
        """ Reads Demux Trace Header """
        return self.read_and_parse_header(start_byte, 20, "_demuxTraceHeader")

        
    def read_trace_header_extension(self, start_byte):
        """ Reads Trace Header Extension """
        return self.read_and_parse_header(start_byte, 32, "_traceHeaderExtension")


    def read_other_header(self, start_byte):
        """ Reads other / optional header extensions given in SEG-D 3.0 documentation """
        if not self.file:
            return None
        if not self.major_version:
            self.get_segd_version()
        if int(self.major_version) < 3:
            return None
    
        trace_header = self.read_32byte_header(start_byte)
        byte32 = self.read_unsigned_binary(trace_header, 32, 32)
        hex_key = hex(byte32)
    
        # Switch to 96-byte read for scan type header and position blocks header
        if hex_key in {"0x30", "0x50"}:
            trace_header = self.read_96byte_header(start_byte)

        # Skip continuations of scan type header and position blocks header
        if hex_key in {"0x31", "0x32", "0x51", "0x52"}:
            if self.verbose:
                print(f"Skipping continuation header: {hex_key}")
            return None
    
        try:
            header_spec = self.read_header_spec_from_json(self.header_specs, hex_key)
        except ValueError:
            if self.verbose:
                print(f"Header type unrecognised: {hex_key}. Skipping.")
            return None
    
        header_fields = self.parse_header_fields_from_spec(trace_header, header_spec)
    
        if self.verbose:
            print(f"Header type: {hex_key}")
            for field, value in header_fields.items():
                print(f"{field}: {value}")
    
        return header_fields



    def read_32bit_IEEE_trace_data(self, start_byte, num_samples):
        """
        Read trace data from a SEG-D file (format code = 8058).
        Trace format used by Sercel and SmartSolo.
        
        Parameters:
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
        

    def read_24bit_trace_data(self, start_byte, num_samples):
        """
        Read 24-bit 2's complement trace data from a SEG-D file (format code = 8036). 
        Trace format used by Stryde rev 3.0.
        
        Parameters:
        - start_byte: The starting byte of the trace data (1-based index).
        - num_samples: Number of samples in the trace.
        
        Returns:
        - A list of decoded trace samples as 24-bit signed integers.
        """        
    
        if not self.file:
            return None
    
        self.file.seek(start_byte)  # Move to the start of trace data
        total_trace_bytes = num_samples * 3  # Total size of trace data in bytes (24 bits per sample)
        relevant_bytes = self.file.read(total_trace_bytes)  # Read the required number of bytes for trace
        
        # Convert each 3-byte sequence to a signed 24-bit integer
        trace_samples = []
        for i in range(num_samples):
            # Read 3 bytes for each sample
            sample_bytes = relevant_bytes[i*3:(i+1)*3]
            
            # Convert 3 bytes to an integer (big-endian)
            int_value = int.from_bytes(sample_bytes, byteorder='big', signed=False)
            
            # Convert to 2's complement signed integer if the value is negative
            if int_value >= 0x800000:  # 0x800000 is the threshold for 24-bit negative numbers
                int_value -= 0x1000000  # Adjust for 24-bit 2's complement
            
            trace_samples.append(int_value)
    
        return trace_samples


    def get_segd_version(self):
        """
        Identifies file SEG-D file version
        """
        if not self.file:
            try:
                self.open_file()
            except:
                print("Can't open file.")

        # Get version info from bytes 11 and 12 of general header 2:
        gen_head2 = self.read_32byte_header(32)
        self.major_version = self.read_unsigned_binary(gen_head2, 11, 11)
        self.minor_version = self.read_unsigned_binary(gen_head2, 12, 12)

        if self.verbose:
            print("SEG-D version number: " + str(self.major_version) + "." + str(self.minor_version))

        return "SEG-D version number: " + str(self.major_version) + "." + str(self.minor_version)


    def clear_segd_version(self):
        self.major_version = None
        self.minor_version = None

    
    def package_has_versioned_header_data(self):
        try:
            contents = pkg_resources.files("SEGD_reader").joinpath("json").iterdir()
            return any(entry.name == f"rev_{self.major_version}_{self.minor_version}" for entry in contents)
        except Exception:
            return False
    

    def get_general_headers_only(self, force_version=None):
        
        if not self.file:
            try:
                self.open_file()
            except:
                print("Can't open file.")

        self.clear_segd_version()

        if force_version:
            try:
                self.major_version, self.minor_version = [int(ver_num) for ver_num in str(force_version).split('.')]
            except:
                print(f"Forced version '{force_version}' not recognised. Reading version from SEG-D file instead.")
        
        if not self.major_version:
            self.get_segd_version()

        
        # Check if the directory exists within 'json'        
        header_spec_available = self.package_has_versioned_header_data()
        if not header_spec_available:
            print(f"Version {self.major_version}.{self.minor_version} not yet supported. May get weird values.")
            if self.major_version < 3:
                # Use rev 2.1 header structure if older version
                self.major_version = 2
                self.minor_version = 1
            else:
                # Use rev 3.0 header structure if newer version
                self.major_version = 3
                self.minor_version = 0

        version_str = f"{self.major_version}_{self.minor_version}"
        json_rel_path = f"json/rev_{version_str}/header_specs_rev_{version_str}.json"
        self.load_header_specs(pkg_resources.files("SEGD_reader").joinpath(json_rel_path))

        out_dict = {}

        if self.verbose:
            print("*** Reading file headers ***")
        # Read first three general headers and extract key info
        out_dict['general_header1'] = self.read_general_header1()
        out_dict['general_header2'] = self.read_general_header2()
        out_dict['general_header3'] = self.read_general_header3() # Doesn't read anything for versions < 3.0
        
        return out_dict
    
    
    def get_data_plus_key_headers(self, force_version=None):
        """
        Stores key headers and trace data in dictionary, that can be used to create obspy stream etc.
        """
        if not self.file:
            try:
                self.open_file()
            except:
                print("Can't open file.")

        if force_version:
            try:
                self.major_version, self.minor_version = [int(ver_num) for ver_num in str(force_version).split('.')]
            except:
                print(f"Forced version '{force_version}' not recognised. Reading version from SEG-D file instead.")
        
        if not self.major_version:
            self.get_segd_version()
        
        # Check if the directory exists within 'json'        
        header_spec_available = self.package_has_versioned_header_data()
        if not header_spec_available:
            print(f"Version {self.major_version}.{self.minor_version} not yet supported. May get weird values.")
            if self.major_version < 3:
                # Use rev 2.1 header structure if older version
                self.major_version = 2
                self.minor_version = 1
            else:
                # Use rev 3.0 header structure if newer version
                self.major_version = 3
                self.minor_version = 0

        version_str = f"{self.major_version}_{self.minor_version}"
        json_rel_path = f"json/rev_{version_str}/header_specs_rev_{version_str}.json"
        self.load_header_specs(pkg_resources.files("SEGD_reader").joinpath(json_rel_path))
        
        out_dict = {}
        out_dict['file_header'] = {}

        if self.verbose:
            print("*** Reading file headers ***")
        # Read first three general headers and extract key info
        gen_head1 = self.read_general_header1()
        gen_head2 = self.read_general_header2()
        gen_head3 = self.read_general_header3() # Doesn't read anything for versions < 3.0

        # Read additional general headers

        fileNumber = gen_head1['fileNumber']
        if fileNumber == "ffff":
            fileNumber = gen_head2['expandedFileNumber']

        duration = gen_head1['recordLength']
        if 'f' in duration:
            duration = gen_head2['extendedRecordLength']
        if int(self.major_version) < 3:
            duration *= 1e3 # Duration given in milliseconds in earlier SEG-D versions
        
        samplingInterval = gen_head1['baseScanInterval']
        if int(self.major_version) < 3:
            samplingInterval *= 1e3 # Given in milliseconds in earlier SEG-D versions
        else:
            if samplingInterval == 255/16:
                samplingInterval = gen_head2['dominantSamplingInterval']
        
        out_dict['file_header']['fileNumber'] = fileNumber
        out_dict['file_header']['SEGDVersion'] = str(gen_head2['majorSEGDrevisionNumber']) + "." + str(gen_head2['minorSEGDrevisionNumber'])
        try:
            out_dict['file_header']['manufacturer'] = MANUFACTURERS_CODE[int(gen_head1['manufacturersCode'])]
        except KeyError:
            out_dict['file_header']['manufacturer'] = str(gen_head1['manufacturersCode']) + " (not defined in SEG-D documentation)"
        out_dict['file_header']['manufacturer_serial'] = str(gen_head1['manunfacturersSerialNumber'])
        out_dict['file_header']['year'] = gen_head1['year']
        out_dict['file_header']['jday'] = gen_head1['day']
        out_dict['file_header']['hour'] = gen_head1['hour']
        out_dict['file_header']['minute'] = gen_head1['minute']
        out_dict['file_header']['second'] = gen_head1['second']
        out_dict['file_header']['duration_microsec'] = duration
        out_dict['file_header']['dominant_sampling_interval_microsec'] = samplingInterval
        
        if int(self.major_version) >= 3:
            timezero = self.segd_timestamp_to_UTCDateTime(gen_head3['timeZero'])
        else:
            timezero = UTCDateTime("20" + gen_head1['year'].zfill(2) + gen_head1['day'].zfill(3) + "T" + gen_head1['hour'].zfill(2) + gen_head1['minute'].zfill(2) + gen_head1['second'].zfill(2))
        
        out_dict['file_header']['record_timezero_utc'] = timezero
        out_dict['file_header']['trace_format_code'] = gen_head1['formatCode']

        # From SEG-D version 3.0
        if int(self.major_version) >= 3:        
            headerSize = gen_head3['headerSize']
            out_dict['file_header']['headerSize_bytes'] = headerSize

            # Loop through remaining general header blocks while i < headerSize
            i = 3*32 # Start after first three general headers:
            num_channel_sets = 0
            while i < headerSize:
                tmp = self.read_other_header(i) # Returns None if unrecognised header, otherwise dict
                if type(tmp) is dict:
                    if tmp['headerBlockType'] == 48:
                        if (tmp['headerBlockType2'] == 49) and (tmp['headerBlockType3'] == 50):
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
                
                #### NEED TO FIGURE OUT WHAT TO DO WHEN 'Number of scan type headers' is a weird value
                # Looks like seismic data always has '1' of these headers though so can just output data at this point
                if int(out_dict['channelSet_' + str(channel_set + 1)]['description']['scanTypeNumber']) > 1:
                    continue
                    
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
                            tmp = self.read_trace_header_extension(start_byte)
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
                                    out_dict['channelSet_' + str(channel_set + 1)]['traceData'][line_num][point_num]['trace_header']['timeZero_utc'] = self.segd_timestamp_to_UTCDateTime(tmp['timeZero'])
                                elif tmp['headerBlockType'] == 65:
                                    out_dict['channelSet_' + str(channel_set + 1)]['traceData'][line_num][point_num]['trace_header']['sensorSensitivity'] = tmp['sensorSensitivity']
                                    out_dict['channelSet_' + str(channel_set + 1)]['traceData'][line_num][point_num]['trace_header']['serialNumber'] = tmp['serialNumber']
                                elif tmp['headerBlockType'] == 64:
                                    out_dict['channelSet_' + str(channel_set + 1)]['traceData'][line_num][point_num]['trace_header']['timeDriftBlock'] = tmp
                                elif tmp['headerBlockType'] == 214:
                                    out_dict['channelSet_' + str(channel_set + 1)]['traceData'][line_num][point_num]['trace_header']['latitude'] = tmp['latitude']
                                    out_dict['channelSet_' + str(channel_set + 1)]['traceData'][line_num][point_num]['trace_header']['longitude'] = tmp['longitude']
                                    out_dict['channelSet_' + str(channel_set + 1)]['traceData'][line_num][point_num]['trace_header']['elevation'] = tmp['elevation']
                        start_byte += 32
                    # Get trace data:
                
                    if out_dict['file_header']['trace_format_code'] == '8058':
                        trace_data = self.read_32bit_IEEE_trace_data(start_byte, out_dict['channelSet_' + str(channel_set + 1)]['traceData'][line_num][point_num]['trace_header']['numSamplesPerTrace'])
                        start_byte += out_dict['channelSet_' + str(channel_set + 1)]['traceData'][line_num][point_num]['trace_header']['numSamplesPerTrace'] * 4
                    elif out_dict['file_header']['trace_format_code'] == '8036':
                        trace_data = self.read_24bit_trace_data(start_byte, out_dict['channelSet_' + str(channel_set + 1)]['traceData'][line_num][point_num]['trace_header']['numSamplesPerTrace'])
                        start_byte += out_dict['channelSet_' + str(channel_set + 1)]['traceData'][line_num][point_num]['trace_header']['numSamplesPerTrace'] * 3
                    else:
                        print('Format code for trace data (' + out_dict['file_header']['trace_format_code'] + ') not yet supported. Exiting.')
                        return None
    
                    # If seismic data, keep, else skip.
                    if out_dict['channelSet_' + str(channel_set + 1)]['description']['channelType'][2:] == '10':
                        out_dict['channelSet_' + str(channel_set + 1)]['traceData'][line_num][point_num]['trace_data'] = trace_data
                            
                    else:
                        continue

        else:
            # For SEG-D version 2.1
            additional_gen_head_blocks = int(gen_head1['generalHeaderBlocks'])       # Number of additional general header blocks, including general header 2
            scan_types_per_record = int(gen_head1['scanTypesPerRecord'])             # Number of scan types per record (loop through these?)
            channel_sets_per_scan_type = int(gen_head1['channelSetsPerScanType'])    # Number of channel sets per scan type (nested loop?)
            num_additional_headers_per_scan_type = int(gen_head1['skewBlocks'])      # Number of 32 byte headers following each scan type header
            extended_headers = gen_head1['extendedHeaderBlocks']                     # Number of 32 byte extended headers (additional equipment info)
            if 'f' in extended_headers:
                extended_headers = int(gen_head2['extendedHeaderBlocks'])
            else:
                extended_headers = int(extended_headers)
            external_headers = gen_head1['externalHeaderBlocks']                     # Number of 32 byte external headers (additional user supplied info)
            if 'f' in external_headers:
                external_headers = int(gen_head2['externalHeaderBlocks'])
            else:
                external_headers = int(external_headers)

            # Skip additional general header blocks and loop through scan type headers:
            start_byte = (1 + additional_gen_head_blocks) * 32
            for _ in range(scan_types_per_record):
                for _ in range(channel_sets_per_scan_type):
                    channel_set_header = self.read_scan_type_header(start_byte)
                    if (int(channel_set_header['scanTypeNumber']) > 0) and (int(channel_set_header['channelSetNumber']) > 0):
                        scan_type_number = channel_set_header['scanTypeNumber']
                        if scan_type_number not in out_dict:
                            out_dict[scan_type_number] = {}
                        channel_set_number = channel_set_header['channelSetNumber']
                        if channel_set_number not in out_dict[scan_type_number]:
                            out_dict[scan_type_number][channel_set_number] = {}
                            out_dict[scan_type_number][channel_set_number]['description'] = channel_set_header
                            out_dict[scan_type_number][channel_set_number]['traceData'] = {}
                            
                    start_byte += 32
            
                # Not sure what to do with skew blocks yet (the number of these is zero in SmartSolo/Sercel example files)
                if num_additional_headers_per_scan_type > 0:
                    for skew_block in range(num_additional_headers_per_scan_type):
                        skew_block_header = self.read_scan_type_header(start_byte) 
                        start_byte += 32

            # Skip past extended and external headers for now as not sure what info is in them (user defined):
            start_byte += (extended_headers + external_headers) * 32

            for dict_key in list(out_dict.keys()):
                if dict_key == 'file_header':
                    continue
                else:
                    for channel_set in list(out_dict[dict_key].keys()):
                        num_traces = int(out_dict[dict_key][channel_set]['description']['numberOfChannelsThisSet'])
                        for trace in range(num_traces):
                            tmp = self.read_demux_trace_header(start_byte)
                            start_byte += 20
                            trace_num = tmp['traceNumber']
                            num_trace_extension_headers = tmp['traceHeaderExtension']
                            for trace_extension_header_num in range(int(num_trace_extension_headers)):
                                # First trace extension header has sensor info, plus line and point numbers
                                if trace_extension_header_num == 0:
                                    trace_header = self.read_trace_header_extension(start_byte)
                                    sensor_type = trace_header['sensorType']
                                    if (sensor_type > 0) and (sensor_type < 10) and (trace_header['numberOfSamplesPerTrace'] > 0):
                                        line_num = trace_header['receiverLineNumber']
                                        if line_num < 0:
                                            line_num = trace_header['extendedReceiverLineNumberInteger']
                                            line_num = 'line_' + str(line_num)
                                        point_num = trace_header['receiverPointNumber']
                                        if point_num < 0:
                                            point_num = trace_header['extendedReceiverPointNumberInteger']
                                            point_num = 'point_' + str(point_num)
                                        if line_num not in out_dict[dict_key][channel_set]['traceData']:
                                            out_dict[dict_key][channel_set]['traceData'][line_num] = {}
                                        if point_num not in out_dict[dict_key][channel_set]['traceData'][line_num]:
                                            out_dict[dict_key][channel_set]['traceData'][line_num][point_num] = {}
                                        if trace_num not in out_dict[dict_key][channel_set]['traceData'][line_num][point_num]:
                                            out_dict[dict_key][channel_set]['traceData'][line_num][point_num][trace_num] = {}
                                        out_dict[dict_key][channel_set]['traceData'][line_num][point_num][trace_num]['trace_header'] = {}
                                        out_dict[dict_key][channel_set]['traceData'][line_num][point_num][trace_num]['trace_header']['sensorType_code'] = sensor_type
                                        out_dict[dict_key][channel_set]['traceData'][line_num][point_num][trace_num]['trace_header']['numSamplesPerTrace'] = trace_header['numberOfSamplesPerTrace']

                                # Remaining trace header extensions are user/manufacturer defined. 
                                
                                # Sercel headers
                                if int(gen_head1['manufacturersCode']) == 13:
                                    
                                    if trace_extension_header_num == 4:
                                        trace_ext_head = self.read_32byte_header(start_byte)
                                        out_dict[dict_key][channel_set]['traceData'][line_num][point_num][trace_num]['trace_header']['longitude'] = self.read_ieee_double_float(trace_ext_head, 9, 16)
                                        out_dict[dict_key][channel_set]['traceData'][line_num][point_num][trace_num]['trace_header']['latitude'] = self.read_ieee_double_float(trace_ext_head, 17, 24)
                                        out_dict[dict_key][channel_set]['traceData'][line_num][point_num][trace_num]['trace_header']['elevation'] = self.read_ieee_float(trace_ext_head, 29, 32)
                                    
                                    if trace_extension_header_num == 5:
                                        trace_ext_head = self.read_32byte_header(start_byte)
                                        out_dict[dict_key][channel_set]['traceData'][line_num][point_num][trace_num]['trace_header']['serialNumber'] = self.read_unsigned_binary(trace_ext_head, 2, 4)
                                        out_dict[dict_key][channel_set]['traceData'][line_num][point_num][trace_num]['trace_header']['sensorSensitivity'] = self.read_ieee_float(trace_ext_head, 21, 24)
                                        

                                # TO DO: add other manufacture headers

                                start_byte += 32
    
                            # Get trace data
                            if out_dict['file_header']['trace_format_code'] == '8058':
                                trace_data = self.read_32bit_IEEE_trace_data(start_byte, out_dict[dict_key][channel_set]['traceData'][line_num][point_num][trace_num]['trace_header']['numSamplesPerTrace'])
                                start_byte += out_dict[dict_key][channel_set]['traceData'][line_num][point_num][trace_num]['trace_header']['numSamplesPerTrace'] * 4
                            elif out_dict['file_header']['trace_format_code'] == '8036':
                                trace_data = self.read_24bit_trace_data(start_byte, out_dict[dict_key][channel_set]['traceData'][line_num][point_num][trace_num]['trace_header']['numSamplesPerTrace'])
                                start_byte += out_dict['channelSet_' + str(channel_set + 1)]['traceData'][line_num][point_num][trace_num]['trace_header']['numSamplesPerTrace'] * 3
                            else:
                                print('Format code for trace data (' + out_dict['file_header']['trace_format_code'] + ') not yet supported. Exiting.')
                                return None
    
                            out_dict[dict_key][channel_set]['traceData'][line_num][point_num][trace_num]['trace_data'] = trace_data


            # Collapse scan type number and channel set number to match levels of SEG-D rev 3.0 out_dict:
            new_out_dict = {}
            
            # Preserve the 'file_header' key without modification
            if 'file_header' in out_dict:
                new_out_dict['file_header'] = out_dict['file_header']
        
            # Iterate over the main keys (except 'file_header') to collapse nested levels
            for first_key, first_value in out_dict.items():
                if first_key == 'file_header':
                    continue  # Skip 'file_header' as it should not be collapsed
        
                for second_key, second_value in first_value.items():
                    # Form the new key by combining the first and second level keys
                    combined_key = f"{first_key}_{second_key}"
                    
                    # Add the entries under the combined key in the collapsed dictionary
                    new_out_dict[combined_key] = second_value

            out_dict = new_out_dict
                
        return out_dict


def SEG_D_to_stream(filelist, convert_to_int = True, use_descale_multiplier = True, apply_sensitivity_correction = False, sensitivity_value = None, serial_to_station_name_dict = None, network_code = 'AA', remove_gaps = False, remove_stations_with_zero_data = False, get_avg_lat_lon_ele = False, reader_verbose = False, forced_segd_version = None, debug=False):
    '''
    Reads a Sercel SEG-D file and returns an obspy stream
    Currently supports SEG-D revisions 2.1 and 3.0
    '''

    # Initiate obspy stream for all data in filelist
    st = Stream()
    
    if type(filelist) is str:
        filelist = [filelist]

    if get_avg_lat_lon_ele:
        lat_lon_ele_dict = {}

    for filepath in filelist:
        # Initiate SEG_D_Reader:
        reader = SEG_D_Reader(filepath, verbose=reader_verbose)
    
        # Open file:
        reader.open_file()
    
        # Get key header info and trace data:
        data = reader.get_data_plus_key_headers(force_version=forced_segd_version)

        if data is None:
            return None
    
        channel_sets = [key for key in list(data.keys()) if 'file_header' not in key]
        channel_sets.sort()
        line_nums = [list(data[chan_set]['traceData'].keys()) for chan_set in channel_sets]
    
        for chan_set, line_num_list in zip(channel_sets, line_nums):
            if len(line_num_list) < 1:
                continue
            for line_num in line_num_list:
                for point_num in list(data[chan_set]['traceData'][line_num].keys()):
                    if 'trace_data' not in data[chan_set]['traceData'][line_num][point_num]:
                        if 'trace_header' in data[chan_set]['traceData'][line_num][point_num]:
                            continue
                        # SEG-D v2.1 code has additional level to dict (trace_num), so check if this exists:
                        if len(list(data[chan_set]['traceData'][line_num][point_num].keys())) > 0:
                            # Extract the first trace section header as the base header for concatenated trace
                            first_trace = next(iter(data[chan_set]['traceData'][line_num][point_num]))
                            # Check if more than one trace for this sensor (point_num):
                            # if len(list(data[chan_set]['traceData'][line_num][point_num].keys())) > 1:
                            # Initialize an empty list to hold trace data for concatenation
                            combined_trace_data = []
                            combined_header = data[chan_set]['traceData'][line_num][point_num][first_trace]['trace_header'].copy()
                            # Initialize the total sample count
                            total_trace_samples = 0
                            # Iterate over each trace_num ('0001', '0002', etc.) to gather data
                            for trace_key, trace_data in data[chan_set]['traceData'][line_num][point_num].items():
                                # Add the trace data to the combined list
                                combined_trace_data.extend(trace_data['trace_data'])
                                # Update the total sample count
                                total_trace_samples += trace_data['trace_header']['numSamplesPerTrace']
                            
                            # Update the combined header's 'numberOfSamplesInTrace' field
                            combined_header['numSamplesPerTrace'] = total_trace_samples
                            # Overwrite the dictionary for this point_num to contain only the combined data
                            data[chan_set]['traceData'][line_num][point_num] = {'trace_header': combined_header, 'trace_data': combined_trace_data}
                        else:
                            continue
                    trace_h, trace_d = data[chan_set]['traceData'][line_num][point_num]['trace_header'], np.array(data[chan_set]['traceData'][line_num][point_num]['trace_data'], dtype=np.float32)

                    if remove_gaps:
                        if 'timeZero_utc' in trace_h:
                            if trace_h['timeZero_utc'] <= UTCDateTime("1980-01-07T00:00:00.000000Z"):
                                continue

                    # if 'serialNumber' not in trace_h:
                    #     trace_h['serialNumber'] = str(line_num) + '_' + str(point_num)
                    
                    tr = Trace(trace_d)

                    if serial_to_station_name_dict is not None:
                        if 'serialNumber' in trace_h:
                            if trace_h['serialNumber'] in serial_to_station_name_dict:
                                if 'network' in serial_to_station_name_dict[trace_h['serialNumber']]:
                                    tr.stats.network = serial_to_station_name_dict[trace_h['serialNumber']]['network_code']
                                else:
                                    tr.stats.network = network_code
                                tr.stats.station = serial_to_station_name_dict[trace_h['serialNumber']]['station_code']
                    else:
                        tr.stats.network = network_code
                        if 'serialNumber' in trace_h:
                            tr.stats.station = str(trace_h['serialNumber'])
                        else:
                            tr.stats.station = str(line_num) + '_' + str(point_num)

                    if 'samplingInterval' in data[chan_set]['description']:
                        sample_rate = data[chan_set]['description']['samplingInterval'] / 1e6 # Given in microseconds
                    else:
                        sample_rate = data['file_header']['dominant_sampling_interval_microsec'] / 1e6 # Given in microseconds
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

                    tr.stats.starttime = data['file_header']['record_timezero_utc']
                    tr.stats._units = "raw"
                    tr.stats.segd = {}
                    if 'serialNumber' in trace_h:
                        tr.stats.segd['serialNumber'] = str(trace_h['serialNumber'])
                    tr.stats.segd.update(data['file_header'])
                    tr.stats.segd.update(data[chan_set]['description'])

                    tr.stats.segd['sensorType_code'] = instrument_code
                    tr.stats.segd['sensor_type'] = SENSOR_CODE[instrument_code]

                    tr.stats._descale_multiplier = 'NOT CORRECTED'
                    tr.stats._sensitivity = 'NOT CORRECTED'

                    if (use_descale_multiplier) & ('descaleMultiplier' in list(tr.stats.segd.keys())):
                        tr_descale_mult = np.float32(float(tr.stats.segd['descaleMultiplier']))
                        if tr_descale_mult != 0:
                            tr.data = tr.data * tr_descale_mult
                            tr.stats._units = 'V'
                            tr.stats._descale_multiplier = 'CORRECTED'
                        del(tr_descale_mult)
                        
                    if 'physicalUnit' in tr.stats.segd:
                        tr.stats.segd['physicalUnit'] = PHYSICAL_UNIT_CODE[int(tr.stats.segd['physicalUnit'])]
                        
                    if 'sensorSensitivity' in trace_h:
                        tr.stats.segd['sensitivity'] = trace_h['sensorSensitivity']
                    elif sensitivity_value is not None:
                        tr.stats.segd['sensitivity'] = sensitivity_value
                    else:
                        tr.stats._sensitivity = 'NOT CORRECTED (SENSITIVITY UNKNOWN)'
                    
                    # Apply sensitivity correction if needed/possible
                    if ('sensitivity' in tr.stats.segd) & (apply_sensitivity_correction) & (use_descale_multiplier) & ('descaleMultiplier' in list(tr.stats.segd.keys())):
                        if tr.stats.segd['sensitivity'] != 0:
                            tr.data = tr.data / tr.stats.segd['sensitivity']
                            tr.stats._sensitivity = 'CORRECTED'
                            if 'physicalUnit' in tr.stats.segd:
                                tr.stats._units = tr.stats.segd['physicalUnit']
                                tr.stats.segd['physicalUnit'] = tr.stats.segd['physicalUnit'] + " (sensitivity correction applied)"
                    else:
                        if 'physicalUnit' in tr.stats.segd:
                            tr.stats.segd['physicalUnit'] = tr.stats.segd['physicalUnit'] + " (sensitivity correction NOT applied)"
                    if 'latitude' in trace_h:
                        tr.stats.segd['latitude'] = trace_h['latitude']
                        if get_avg_lat_lon_ele:
                            if tr.stats.station not in list(lat_lon_ele_dict.keys()):
                                lat_lon_ele_dict[tr.stats.station] = {}
                            if 'latitude' not in list(lat_lon_ele_dict[tr.stats.station].keys()):
                                lat_lon_ele_dict[tr.stats.station]['latitude'] = []
                            lat_lon_ele_dict[tr.stats.station]['latitude'].append(np.float32(trace_h['latitude']))
                    if 'longitude' in trace_h:
                        tr.stats.segd['longitude'] = trace_h['longitude']
                        if get_avg_lat_lon_ele:
                            if tr.stats.station not in list(lat_lon_ele_dict.keys()):
                                lat_lon_ele_dict[tr.stats.station] = {}
                            if 'longitude' not in list(lat_lon_ele_dict[tr.stats.station].keys()):
                                lat_lon_ele_dict[tr.stats.station]['longitude'] = []
                            lat_lon_ele_dict[tr.stats.station]['longitude'].append(np.float32(trace_h['longitude']))
                    if 'elevation' in trace_h:
                        tr.stats.segd['elevation'] = trace_h['elevation']
                        if get_avg_lat_lon_ele:
                            if tr.stats.station not in list(lat_lon_ele_dict.keys()):
                                lat_lon_ele_dict[tr.stats.station] = {}
                            if 'elevation' not in list(lat_lon_ele_dict[tr.stats.station].keys()):
                                lat_lon_ele_dict[tr.stats.station]['elevation'] = []
                            lat_lon_ele_dict[tr.stats.station]['elevation'].append(np.float32(trace_h['elevation']))
                    # st.append(tr)
                    st += tr
                    
        reader.close_file()

    # Remove any stations with all zero data:
    if remove_stations_with_zero_data:
        stations = list({tr.stats.station for tr in st})
        for sta in stations:
            st_check = st.select(station=sta)
            if np.all([np.all((tr.data == 0) | (np.isnan(tr.data))) for tr in st_check]):
                for tr in st_check:
                    st.remove(tr)
                    
    # Check if all traces can be converted to int
    convert_to_int = convert_to_int and np.all([(np.mod(tr.data, 1) == 0) for tr in st])
    
    if convert_to_int:
        for tr in st:
            tr.data = tr.data.astype(np.int32)
    
    try:
        st.merge() # Merge consecutive files
    except Exception as e:
        # Handle the exception or ignore it
        print(f"Ignoring error: {e}")

    # Assign average lat, lon, ele:
    if get_avg_lat_lon_ele:
        for tr in st:
            sta = tr.stats.station
            if 'loc' not in list(tr.stats.keys()):
                tr.stats.loc = {}
            tr.stats.loc['latitude'] = np.mean(lat_lon_ele_dict[tr.stats.station]['latitude'])
            tr.stats.loc['longitude'] = np.mean(lat_lon_ele_dict[tr.stats.station]['longitude'])
            tr.stats.loc['elevation'] = np.mean(lat_lon_ele_dict[tr.stats.station]['elevation'])

    # Make sure 3C nodes have all three components (fill with zeros if not? Or masked array?)
    if remove_gaps:
        if serial_to_station_name_dict is not None:
            stations = list({tr.stats.station for tr in st})
            for sta in stations:
                st_check = st.copy().select(station=sta)
                station_components = defaultdict(list)
                for node_serial, details in serial_to_station_name_dict.items():
                    station_code = details['station_code']
                    component = details['component']
                    station_components[station_code].append({'component': component, 'node_serial': node_serial})
                for station_node in list(station_components[sta]):
                    if np.all([(tr.stats.component != station_node['component']) for tr in st_check]):
                        new_trace = st_check.copy()[0] # Copy first component from station
                        new_trace.data = np.zeros(new_trace.data.shape) # Replace data with all zeroes
                        new_trace.stats.channel = new_trace.stats.channel[:2] + station_node['component'] # Update component name
                        new_trace.stats.segd.serialNumber = station_node['node_serial'] # Update node serial number
                        # Add new_trace on to st_check and st
                        st_check += new_trace.copy()
                        st += new_trace.copy()

            st.sort()

    if get_avg_lat_lon_ele:
        return st, lat_lon_ele_dict
        
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

MANUFACTURERS_CODE = {
    1: 'Alpine Geophysical Associates, Inc. (Obsolete)',
    2: 'Applied Magnetics Corporation (See 09)',
    3: 'Western Geophysical Exploration Products (formerly Litton Resources Systems)',
    4: 'SIE, Inc. (Obsolete)',
    5: 'Dyna-Tronics Mfg. Corporation (Obsolete)',
    6: 'Electronic Instrumentation, Inc. (Obsolete)',
    7: 'Halliburton Geophysical Services, Inc. (formerly, Electro-Technical Labs, Div. of Geosource, Inc.)',
    8: 'Fortune Electronics, Inc. (Obsolete)',
    9: 'Geo Space Corporation',
    10: 'Leach Corporation (Obsolete)',
    11: 'Metrix Instrument Co. (Obsolete)',
    12: 'Redcor Corporation (Obsolete)',
    13: "Sercel (Societe d'Etudes, Recherches Et Constructions Electroniques)",
    14: 'Scientific Data Systems (SDS), (Obsolete)',
    15: 'Texas Instruments, Inc.',
    17: 'GUS Manufacturing, Inc.',
    18: 'Input/Output, Inc.',
    19: 'Geco-Prakla',
    20: 'Fairfield Industries, Incorporated',
    22: 'Geco-Prakla',
    31: 'Japex Geoscience Institute',
    32: 'Halliburton Geophysical Services, Inc.',
    33: 'Compuseis, Inc.',
    34: 'Syntron, Inc.',
    35: 'Syntron Europe Ltd.',
    36: 'Opseis',
    39: 'Grant Geophysical',
    40: 'Geo-X',
    41: 'PGS Inc.',
    42: 'Seamap UK Ltd.',
    43: 'Hydroscience',
    44: 'JSC',
    45: 'Fugro',
    46: 'ProFocus Systems AS',
    47: 'Optoplan AS',
    48: 'Wireless Seismic Inc.',
    49: 'AutoSeis',
    50: 'INOVA Geophysical, Inc.',
    51: 'Verif-i Ltd.',
    52: 'Troika International',
    53: 'MagSeis AS',
    54: 'Seismic Instruments, Inc.',
    55: 'TGS',
    56: 'Hewlett Packard Co',
    57: 'Modern Seismic Technology, LLC'
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