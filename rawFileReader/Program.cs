// --------------------------------------------------------------------------------------------------------------------
// <copyright file="Program.cs" company="Thermo Fisher Scientific">
//   Copyright © Thermo Fisher Scientific. All rights reserved.
// </copyright>
// This source file contains the code for the console application that demonstrates who to
// read our RAW files using the RAWFileReader .Net library.  This example uses only a fraction
// of the functions available in the RawFileReader library.  Please consult the RawFileReader
// documentation for a complete list of methods available.
// 
// This code has been compiled and tested using Visual Studio 2013 Update 5, Visual Studio 
// 2015 Update 3, Visual Studio 2017 Update 2, Visual Studio MAC, and MonoDevelop.  Additional 
// tools used include Resharper 2017.1.  This application has been tested with .Net Framework 4.5.1, 
// 4.5.2, 4.6, 4.6.1, 4.6.2, and 4.7.
//
// Questions about the program can be directed to jim.shofstahl@thermofisher.com
// --------------------------------------------------------------------------------------------------------------------
using Newtonsoft.Json;

namespace RawFileReaderDotNetExample
{
    // These are the libraries necessary to read the Thermo RAW files.  The Interfaces
    // library contains the extension for accessing the scan averaging and background
    // subtraction methods.
    using System;
    using System.Collections.Generic;
    using System.Diagnostics;
    using System.IO;
    using System.Linq;
    using System.Runtime.ExceptionServices;

    using ThermoFisher.CommonCore.Data;
    using ThermoFisher.CommonCore.Data.Business;
    using ThermoFisher.CommonCore.Data.FilterEnums;
    using ThermoFisher.CommonCore.Data.Interfaces;
    using ThermoFisher.CommonCore.MassPrecisionEstimator;
    using ThermoFisher.CommonCore.RawFileReader;

    using Newtonsoft.Json;

    /// <summary>
    /// The object used to store the inclusion/exclusion list items.
    /// </summary>
    public class InclusionListItem : IComparable<InclusionListItem>
    {
        /// <summary>
        /// Gets or sets the mass value for the inclusion/exclusion item.
        /// </summary>
        public double Mass { get; set; }

        /// <summary>
        /// Gets or sets the acquisition parameter for the inclusion/exclusion item.  This
        /// could be an intensity threshold or an isolation window based upon the instrument.
        /// </summary>
        public double Threshold { get; set; }

        /// <summary>
        /// Gets or sets the scan number for the inclusion/exclusion item.
        /// </summary>
        public int ScanNumber { get; set; }

        /// <summary>
        /// Gets or sets the descriptor for the inclusion/exclusion item.
        /// </summary>
        public string Descriptor { get; set; }

        /// <summary>
        /// Gets or sets the if the item is an exclusion item.
        /// </summary>
        public bool IsExclusionItem { get; set; }


        /// <summary>
        /// Initializes a new instance of the <see cref="InclusionListItem"/> class.
        /// </summary>
        public InclusionListItem()
        {
            Mass = 0.0;
            ScanNumber = 0;
            Threshold = 0.0;
            Descriptor = string.Empty;
            IsExclusionItem = false;
        }

        /// <summary>
        /// The compare to method for the inclusion item class.
        /// </summary>
        /// <param name="other">
        /// The other InclusionItem object to compare against this object.
        /// </param>
        /// <returns>
        /// The <see cref="int"/>.
        /// </returns>
        public int CompareTo(InclusionListItem other)
        {
            // Compare the descriptors
            int result = string.CompareOrdinal(Descriptor, other.Descriptor);
            if (result != 0)
            {
                return result;
            }

            // Compare the mass values
            result = Mass.CompareTo(other.Mass);
            if (result != 0)
            {
                return result;
            }

            // Compare the threshold values
            result = Threshold.CompareTo(other.Threshold);
            if (result != 0)
            {
                return result;
            }

            // Compare the scan number values
            result = ScanNumber.CompareTo(other.ScanNumber);
            if (result != 0)
            {
                return result;
            }

            return result;
        }
    }

    /// </summary>
    internal static class Program
    {
        /// <summary>
        /// The main routine for this example program.  This example program only shows how to use the RawFileReader library
        /// in a single-threaded application but our documentation for RawFileReader describes how to use it in a multi-threaded
        /// application.
        /// </summary>
        /// <param name="args">The command line arguments for this program.  The RAW file name should be passed as the first argument</param>
        private static void Main(string[] args)
        {
            // This local variables controls if certain operations are performed. Change any of them to true to read and output that
            // information section from the RAW file.
            bool averageScans = false;
            bool getChromatogram = false;
            bool getInclusionExclusionList = false;
            bool getStatusLog = false;
            bool getTrailerExtra = false;
            bool readAllScans = false;
            bool readScanInformation = false;
            bool getSpectrum = true;

            DirectoryInfo di = new DirectoryInfo("/Users/jameshoward/Documents/Programming/PeakPicker/DesiSequencer/data");
            FileInfo[] files = di.GetFiles("*.RAW");
            foreach (FileInfo file in files)
            {
                Console.WriteLine(file.FullName);

                string filename = file.FullName;
                Console.WriteLine(filename);


                // Get the memory used at the beginning of processing
                Process processBefore = Process.GetCurrentProcess();
                long memoryBefore = processBefore.PrivateMemorySize64 / 1024;

                try
                {
                    // Check to see if the RAW file name was supplied as an argument to the program
                    //string filename = file.FullName;

                    if (string.IsNullOrEmpty(filename))
                    {
                        Console.WriteLine("No RAW file specified!");

                        return;
                    }

                    // Check to see if the specified RAW file exists
                    if (!File.Exists(filename))
                    {
                        Console.WriteLine(@"The file doesn't exist in the specified location - " + filename);

                        return;
                    }

                    // Create the IRawDataPlus object for accessing the RAW file
                    var rawFile = RawFileReaderAdapter.FileFactory(filename);

                    if (!rawFile.IsOpen || rawFile.IsError)
                    {
                        Console.WriteLine("Unable to access the RAW file using the RawFileReader class!");

                        return;
                    }

                    // Check for any errors in the RAW file
                    if (rawFile.IsError)
                    {
                        Console.WriteLine("Error opening ({0}) - {1}", rawFile.FileError, filename);

                        return;
                    }

                    // Check if the RAW file is being acquired
                    if (rawFile.InAcquisition)
                    {
                        Console.WriteLine("RAW file still being acquired - " + filename);

                        return;
                    }

                    // Get the number of instruments (controllers) present in the RAW file and set the 
                    // selected instrument to the MS instrument, first instance of it
                    Console.WriteLine("The RAW file has data from {0} instruments" + rawFile.InstrumentCount);

                    rawFile.SelectInstrument(Device.MS, 1);

                    // Get the first and last scan from the RAW file
                    int firstScanNumber = rawFile.RunHeaderEx.FirstSpectrum;
                    int lastScanNumber = rawFile.RunHeaderEx.LastSpectrum;


                    // Get the start and end time from the RAW file
                    double startTime = rawFile.RunHeaderEx.StartTime;
                    double endTime = rawFile.RunHeaderEx.EndTime;

                    // Get some information from the header portions of the RAW file and display that information.
                    // The information is general information pertaining to the RAW file.
                    Console.WriteLine("General File Information:");
                    Console.WriteLine("   RAW file: " + rawFile.FileName);
                    Console.WriteLine("   Mass resolution: {0:F3} ", rawFile.RunHeaderEx.MassResolution);
                    Console.WriteLine("   Number of scans: {0}", rawFile.RunHeaderEx.SpectraCount);
                    Console.WriteLine("   Scan range: {0} - {1}", firstScanNumber, lastScanNumber);
                    Console.WriteLine("   Time range: {0:F2} - {1:F2}", startTime, endTime);
                    Console.WriteLine("   Mass range: {0:F4} - {1:F4}", rawFile.RunHeaderEx.LowMass, rawFile.RunHeaderEx.HighMass);
                    Console.WriteLine();


                    // Read the first instrument method (most likely for the MS portion of the instrument).
                    // NOTE: This method reads the instrument methods from the RAW file but the underlying code
                    // uses some Microsoft code that hasn't been ported to Linux or MacOS.  Therefore this
                    // method won't work on those platforms therefore the check for Windows.
                    if (Environment.OSVersion.ToString().Contains("Windows"))
                    {
                        var deviceNames = rawFile.GetAllInstrumentNamesFromInstrumentMethod();

                        foreach (var device in deviceNames)
                        {
                            Console.WriteLine("Instrument method: " + device);
                        }

                        Console.WriteLine();
                    }

                    // Display all of the trailer extra data fields present in the RAW file
                    if (getTrailerExtra)
                    {
                        ListTrailerExtraFields(rawFile);
                    }


                    // Get the inclusion/exclusion list
                    if (getInclusionExclusionList)
                    {
                        var inclusionList = GetInclusionExclusionList(rawFile, 1e-5);

                        // Output the saved inclusion/exclusion list
                        int count = 0;

                        foreach (var item in inclusionList)
                        {
                            Console.WriteLine("  {0} - {1}, {2:F4}, {3:F0}, {4}", ++count, item.Descriptor, item.Mass, item.Threshold, item.ScanNumber);
                        }

                        Console.WriteLine();
                    }

                    // Get the number of filters present in the RAW file
                    int numberFilters = rawFile.GetFilters().Count;

                    // Get the scan filter for the first and last spectrum in the RAW file
                    var firstFilter = rawFile.GetFilterForScanNumber(firstScanNumber);
                    var lastFilter = rawFile.GetFilterForScanNumber(lastScanNumber);

                    Console.WriteLine("Filter Information:");
                    Console.WriteLine("   Scan filter (first scan): " + firstFilter.ToString());
                    Console.WriteLine("   Scan filter (last scan): " + lastFilter.ToString());
                    Console.WriteLine("   Total number of filters:" + numberFilters);
                    Console.WriteLine();

                    // Get the BasePeak chromatogram for the MS data
                    if (getChromatogram)
                    {
                        GetChromatogram(rawFile, firstScanNumber, lastScanNumber, true);
                    }

                    // Read the scan information for each scan in the RAW file
                    if (readScanInformation)
                    {
                        ReadScanInformation(rawFile, firstScanNumber, lastScanNumber, true);
                    }

                    // Get a spectrum from the RAW file.  
                    if (getSpectrum)
                    {
                        string filepath = rawFile.FileName;
                        string filePathWithoutExt = Path.ChangeExtension(filepath, null);
                        var FilePath = String.Format("{0}.txt", filePathWithoutExt);
                        //File.AppendAllText(FilePath, "[");
                        for (int i = firstScanNumber; i < lastScanNumber; i++)
                        {
                            
                            GetSpectrum(rawFile, i, true);
                            File.AppendAllText(FilePath, "|");

                        }

                        

                    }

                    // Get a average spectrum from the RAW file for the first 15 scans in the file.  
                    if (averageScans)
                    {
                        GetAverageSpectrum(rawFile, firstScanNumber, lastScanNumber, true);
                    }

                    // Read each spectrum
                    if (readAllScans)
                    {
                        ReadAllSpectra(rawFile, firstScanNumber, lastScanNumber, true);
                    }


                    // Close (dispose) the RAW file
                    Console.WriteLine();
                    Console.WriteLine("Closing " + filename);

                    rawFile.Dispose();
                }
                catch (Exception ex)
                {
                    Console.WriteLine("Error accessing RAWFileReader library! - " + ex.Message);
                }

                // Get the memory used at the end of processing
                Process processAfter = Process.GetCurrentProcess();
                long memoryAfter = processAfter.PrivateMemorySize64 / 1024;

                Console.WriteLine();
                Console.WriteLine("Memory Usage:");
                Console.WriteLine("   Before {0} kb, After {1} kb, Extra {2} kb", memoryBefore, memoryAfter, memoryAfter - memoryBefore);
            }
        }

        /// <summary>
        /// Gets the average spectrum from the RAW file.  
        /// </summary>
        /// <param name="rawFile">
        /// The RAW file being read
        /// </param>
        /// <param name="firstScanNumber">
        /// The first scan to consider for the averaged spectrum
        /// </param>
        /// <param name="lastScanNumber">
        /// The last scan to consider for the averaged spectrum
        /// </param>
        /// <param name="outputData">
        /// The output data flag.
        /// </param>
        private static void GetAverageSpectrum(IRawDataPlus rawFile, int firstScanNumber, int lastScanNumber, bool outputData)
        {
            // Create the mass options object that will be used when averaging the scans

            //MAYBE USE FOR SPECIFIC SCANS


            var options = rawFile.DefaultMassOptions();

            options.ToleranceUnits = ToleranceUnits.ppm;
            options.Tolerance = 5.0;

            // Get the scan filter for the first scan.  This scan filter will be used to located
            // scans within the given scan range of the same type
            var scanFilter = rawFile.GetFilterForScanNumber(firstScanNumber);

            // Get the average mass spectrum for the provided scan range. In addition to getting the
            // average scan using a scan range, the library also provides a similar method that takes
            // a time range.
            var averageScan = rawFile.AverageScansInScanRange(firstScanNumber, lastScanNumber, scanFilter, options);

            if (averageScan.HasCentroidStream)
            {
                Console.WriteLine("Average spectrum ({0} points)", averageScan.CentroidScan.Length);

                // Print the spectral data (mass, intensity values)
                if (outputData)
                {
                    for (int i = 0; i < averageScan.CentroidScan.Length; i++)
                    {
                        Console.WriteLine("  {0:F4} {1:F0}", averageScan.CentroidScan.Masses[0], averageScan.CentroidScan.Intensities[i]);
                    }
                }
            }

            // This example uses a different method to get the same average spectrum that was calculated in the
            // previous portion of this method.  Instead of passing the start and end scan, a list of scans will
            // be passed to the GetAveragedMassSpectrum function.
            List<int> scans = new List<int>(new[] { 1, 6, 7, 9, 11, 12, 14 });

            averageScan = rawFile.AverageScans(scans, options);

            if (averageScan.HasCentroidStream)
            {
                Console.WriteLine("Average spectrum ({0} points)", averageScan.CentroidScan.Length);

                // Print the spectral data (mass, intensity values)
                if (outputData)
                {
                    for (int i = 0; i < averageScan.CentroidScan.Length; i++)
                    {
                        Console.WriteLine("  {0:F4} {1:F0}", averageScan.CentroidScan.Masses[0], averageScan.CentroidScan.Intensities[i]);
                    }
                }
            }

            Console.WriteLine();
        }

        /// <summary>
        /// Reads the base peak chromatogram for the RAW file
        /// </summary>
        /// <param name="rawFile">
        /// The RAW file being read
        /// </param>
        /// <param name="startScan">
        /// Start scan for the chromatogram
        /// </param>
        /// <param name="endScan">
        /// End scan for the chromatogram
        /// </param>
        /// <param name="outputData">
        /// The output data flag.
        /// </param>
        private static void GetChromatogram(IRawDataPlus rawFile, int startScan, int endScan, bool outputData)
        {
            // Define the settings for getting the Base Peak chromatogram
            ChromatogramTraceSettings settings = new ChromatogramTraceSettings(TraceType.BasePeak);

            // Get the chromatogram from the RAW file. 
            var data = rawFile.GetChromatogramData(new IChromatogramSettings[] { settings }, startScan, endScan);

            // Split the data into the chromatograms
            var trace = ChromatogramSignal.FromChromatogramData(data);

            if (trace[0].Length > 0)
            {
                // Print the chromatogram data (time, intensity values)
                Console.WriteLine("Base Peak chromatogram ({0} points)", trace[0].Length);

                if (outputData)
                {
                    for (int i = 0; i < trace[0].Length; i++)
                    {
                        Console.WriteLine("  {0} - {1:F3}, {2:F0}", i, trace[0].Times[i], trace[0].Intensities[i]);
                    }
                }
            }

            Console.WriteLine();
        }

        /// <summary>
        /// Reads the inclusion/exclusion list from the mass spectrometer method in the RAW file
        /// </summary>
        /// <param name="rawFile">
        /// The RAW file being read
        /// </param>
        /// <param name="massTolerance">
        /// The mass tolerance used when associating an inclusion list item with a spectrum
        /// </param>
        /// <returns>
        /// The <see cref="List"/>.
        /// </returns>
        private static List<InclusionListItem> GetInclusionExclusionList(IRawDataPlus rawFile, double massTolerance)
        {
            // Select the MS instrument
            rawFile.SelectInstrument(Device.MS, 1);

            // Get the instrument method item(s) and look for the inclusion/exclusion list
            // which will be flagged by the "Mass List Table"
            List<string> inclusionStrings = new List<string>();

            for (int i = 0; i < rawFile.InstrumentMethodsCount; i++)
            {
                var methodText = rawFile.GetInstrumentMethod(i);

                if (methodText.Contains("Mass List Table"))
                {
                    bool saveLine = false;
                    var splitMethod = methodText.Split(new[] { "\n" }, StringSplitOptions.None);

                    foreach (var line in splitMethod)
                    {
                        if (line.Contains("Mass List Table"))
                        {
                            saveLine = true;
                        }

                        if (line.Contains("End Mass List Table"))
                        {
                            saveLine = false;

                            continue;
                        }

                        if (saveLine && !line.Contains("Mass List Table"))
                        {
                            inclusionStrings.Add(line);
                        }
                    }
                }
            }

            // Create the inclusion/exclusion list 
            List<InclusionListItem> inclusionList = new List<InclusionListItem>();

            // Convert each line from the inclusion/exclusion mass table into InclusionListItem objects
            // and add them to the inclusion/exclusion list.
            foreach (var line in inclusionStrings)
            {
                // Skip the title line
                if (line.Contains("CompoundName"))
                {
                    continue;
                }

                // Split the line into its separate fields
                var fields = line.Split(new[] { "|" }, StringSplitOptions.None);

                if (fields.Length == 4)
                {
                    InclusionListItem inclusionItem = new InclusionListItem() { Descriptor = fields[0] };

                    inclusionItem.Mass = Convert.ToDouble(fields[1]);
                    inclusionItem.Threshold = Convert.ToDouble(fields[2]);

                    inclusionList.Add(inclusionItem);
                }
            }

            // Get the actual scan number for each mass in the inclusion list
            for (int scan = rawFile.RunHeaderEx.FirstSpectrum; scan <= rawFile.RunHeaderEx.LastSpectrum; scan++)
            {
                // Get the scan filter and event for this scan number
                var scanFilter = rawFile.GetFilterForScanNumber(scan);
                var scanEvent = rawFile.GetScanEventForScanNumber(scan);

                // Only consider MS2 scans when looking for the spectrum corr3e
                if (scanFilter.MSOrder == MSOrderType.Ms2)
                {
                    // Get the reaction information in order to get the precursor mass for this spectrum
                    var reaction = scanEvent.GetReaction(0);

                    double precursorMass = reaction.PrecursorMass;
                    double tolerance = precursorMass * massTolerance;

                    // Find the inclusion list item matching this precursor mass
                    foreach (var item in inclusionList)
                    {
                        if (item.Mass >= (precursorMass - tolerance) && (item.Mass <= (precursorMass + tolerance)))
                        {
                            item.ScanNumber = scan;
                            break;
                        }
                    }
                }
            }

            return inclusionList;
        }

        /// <summary>
        /// Gets the spectrum from the RAW file.
        /// </summary>
        /// <param name="rawFile">
        /// The RAW file being read
        /// </param>
        /// <param name="scanNumber">
        /// The scan number being read
        /// </param>
        /// <param name="scanFilter">
        /// The scan filter for that scan
        /// </param>
        /// <param name="outputData">
        /// The output data flag.
        /// </param>
        private static void GetSpectrum(IRawDataPlus rawFile, int scanNumber, bool outputData)
        {
            // Get the scan statistics from the RAW file for this scan number
            var scanStatistics = rawFile.GetScanStatsForScanNumber(scanNumber);

            // Check to see if the scan has centroid data or profile data.  Depending upon the
            // type of data, different methods will be used to read the data.  While the ReadAllSpectra
            // method demonstrates reading the data using the Scan.FromFile method, generating the
            // Scan object takes more time and memory to do, so that method isn't optimum.
            if (scanStatistics.IsCentroidScan && scanStatistics.SpectrumPacketType == SpectrumPacketType.FtCentroid)
            {
                // Get the centroid (label) data from the RAW file for this scan
                var centroidStream = rawFile.GetCentroidStream(scanNumber, false);

                var DataForFile = "";

                // Create dictionary to add masses and intensities
                Dictionary<double, double> DESIMasses = new Dictionary<double, double>();

                // Console.WriteLine("Spectrum (centroid/label) {0} - {1} points", scanNumber, centroidStream.Length);

                // Print the spectral data (mass, intensity, charge values).  Not all of the information in the high resolution centroid 
                // (label data) object is reported in this example.  Please check the documentation for more information about what is
                // available in high resolution centroid (label) data.
                if (outputData)
                {
                    for (int i = 0; i < centroidStream.Length; i++)
                    {
                        //DESIMasses.Add(centroidStream.Masses[i], centroidStream.Intensities[i]);

                        //Console.WriteLine("{0} - {1:F4}, {2:F0}", i, centroidStream.Masses[i], centroidStream.Intensities[i]);
                        //var Data = String.Format("{0}, {1:F4}, {2:F0};", i, centroidStream.Masses[i], centroidStream.Intensities[i]);
                        var Data = String.Format("{0:F4}, {1:F0};", centroidStream.Masses[i], centroidStream.Intensities[i]);
                        DataForFile += Data;
                    }
                }
                Console.WriteLine(DataForFile);

                //string DESIJsonMasses = JsonConvert.SerializeObject(DESIMasses);

                string filepath = rawFile.FileName;
                string filePathWithoutExt = Path.ChangeExtension(filepath, null);



                //var FilePath = String.Format("{0}{1}.json", filePathWithoutExt, scanNumber);
                //var WhichScan = String.Format("Scan {0}", scanNumber);
                var FilePath = String.Format("{0}.txt", filePathWithoutExt);
                //File.AppendAllText(FilePath, WhichScan);
                //File.AppendAllText(FilePath, DESIJsonMasses);
                File.AppendAllText(FilePath, DataForFile);

                Console.WriteLine();
            }
            else
            {
                // Get the segmented (low res and profile) scan data
                var segmentedScan = rawFile.GetSegmentedScanFromScanNumber(scanNumber, scanStatistics);

                Console.WriteLine("Spectrum (normal data) {0} - {1} points", scanNumber, segmentedScan.Positions.Length);
                var DataForFile = "";
               

                // Print the spectral data (mass, intensity values)
                if (outputData)
                {
                    for (int i = 0; i < segmentedScan.Positions.Length; i++)
                    {

                        var Data = String.Format("  {0} - {1:F4}, {2:F0}", i, segmentedScan.Positions[i], segmentedScan.Intensities[i]);
                        // Console.WriteLine(Data);
                        DataForFile += Data;





                    }
                }
               
            }
        }

        /// <summary>
        /// Reads and reports the trailer extra data fields present in the RAW file.
        /// </summary>
        /// <param name="rawFile">
        /// The RAW file
        /// </param>
        private static void ListTrailerExtraFields(IRawDataPlus rawFile)
        {
            // Get the Trailer Extra data fields present in the RAW file
            var trailerFields = rawFile.GetTrailerExtraHeaderInformation();

            // Display each value
            int i = 0;

            Console.WriteLine("Trailer Extra Data Information:");

            foreach (var field in trailerFields)
            {
                Console.WriteLine("   Field {0} = {1} storing data of type {2}", i, field.Label, field.DataType);

                i++;
            }

            Console.WriteLine();
        }


        /// <summary>
        /// Read all spectra in the RAW file.
        /// </summary>
        /// <param name="rawFile">
        /// The raw file.
        /// </param>
        /// <param name="firstScanNumber">
        /// The first scan number.
        /// </param>
        /// <param name="lastScanNumber">
        /// The last scan number.
        /// </param>
        /// <param name="outputData">
        /// The output data flag.
        /// </param>
        [HandleProcessCorruptedStateExceptions]
        private static void ReadAllSpectra(IRawDataPlus rawFile, int firstScanNumber, int lastScanNumber, bool outputData)
        {
            for (int scanNumber = firstScanNumber; scanNumber <= lastScanNumber; scanNumber++)
            {
                try
                {
                    // Get the scan filter for the spectrum
                    var scanFilter = rawFile.GetFilterForScanNumber(scanNumber);

                    if (string.IsNullOrEmpty(scanFilter.ToString()))
                    {
                        continue;
                    }

                    // Get the scan from the RAW file.  This method uses the Scan.FromFile method which returns a
                    // Scan object that contains both the segmented and centroid (label) data from an FTMS scan
                    // or just the segmented data in non-FTMS scans.  The GetSpectrum method demonstrates an
                    // alternative method for reading scans.
                    var scan = Scan.FromFile(rawFile, scanNumber);

                    // If that scan contains FTMS data then Centroid stream will be populated so check to see if it is present.
                    int labelSize = 0;

                    if (scan.HasCentroidStream)
                    {
                        labelSize = scan.CentroidScan.Length;
                    }

                    // For non-FTMS data, the preferred data will be populated
                    int dataSize = scan.PreferredMasses.Length;

                    if (outputData)
                    {
                        Console.WriteLine("Spectrum {0} - {1}: normal {2}, label {3} points", scanNumber, scanFilter.ToString(), dataSize, labelSize);
                    }
                }
                catch (Exception ex)
                {
                    Console.WriteLine("Error reading spectrum {0} - {1}", scanNumber, ex.Message);
                }
            }
        }

        /// <summary>
        /// Reads the general scan information for each scan in the RAW file using the scan filter object and also the
        /// trailer extra data section for that same scan.
        /// </summary>
        /// <param name="rawFile">
        /// The RAW file being read
        /// </param>
        /// <param name="firstScanNumber">
        /// The first scan in the RAW file
        /// </param>
        /// <param name="lastScanNumber">
        /// the last scan in the RAW file
        /// </param>
        /// <param name="outputData">
        /// The output data flag.
        /// </param>
        private static void ReadScanInformation(IRawDataPlus rawFile, int firstScanNumber, int lastScanNumber, bool outputData)
        {
            // Read each scan in the RAW File
            for (int scan = firstScanNumber; scan <= lastScanNumber; scan++)
            {
                // Get the retention time for this scan number.  This is one of two comparable functions that will
                // convert between retention time and scan number.
                double time = rawFile.RetentionTimeFromScanNumber(scan);

                // Get the scan filter for this scan number
                // NOTE: A scan filter can also be created from the filter string using the GetFilterFromString in the
                // IRawDataPlus or IFilterParser class.
                var scanFilter = rawFile.GetFilterForScanNumber(scan);

                // Get the scan event for this scan number
                var scanEvent = rawFile.GetScanEventForScanNumber(scan);

                // Get the ionizationMode, MS2 precursor mass, collision energy, and isolation width for each scan
                if (scanFilter.MSOrder == MSOrderType.Ms2)
                {
                    // Get the reaction information for the first precursor
                    var reaction = scanEvent.GetReaction(0);

                    double precursorMass = reaction.PrecursorMass;
                    double collisionEnergy = reaction.CollisionEnergy;
                    double isolationWidth = reaction.IsolationWidth;
                    double monoisotopicMass = 0.0;
                    int masterScan = 0;
                    var ionizationMode = scanFilter.IonizationMode;
                    var order = scanFilter.MSOrder;

                    // Get the trailer extra data for this scan and then look for the monoisotopic m/z value in the 
                    // trailer extra data list
                    var trailerData = rawFile.GetTrailerExtraInformation(scan);

                    for (int i = 0; i < trailerData.Length; i++)
                    {
                        if (trailerData.Labels[i] == "Monoisotopic M/Z:")
                        {
                            monoisotopicMass = Convert.ToDouble(trailerData.Values[i]);
                        }

                        if ((trailerData.Labels[i] == "Master Scan Number:") || (trailerData.Labels[i] == "Master Scan Number") || (trailerData.Labels[i] == "Master Index:"))
                        {
                            masterScan = Convert.ToInt32(trailerData.Values[i]);
                        }
                    }

                    if (outputData)
                    {
                        Console.WriteLine(
                            "Scan number {0} @ time {1:F2} - Master scan = {2}, Ionization mode={3}, MS Order={4}, Precursor mass={5:F4}, Monoisotopic Mass = {6:F4}, Collision energy={7:F2}, Isolation width={8:F2}",
                            scan, time, masterScan, ionizationMode, order, precursorMass, monoisotopicMass, collisionEnergy, isolationWidth);
                    }
                }
                else if (scanFilter.MSOrder == MSOrderType.Ms)
                {
                    var scanDependents = rawFile.GetScanDependents(scan, 5);

                    if (scanDependents != null)
                    {
                        Console.WriteLine(
                            "Scan number {0} @ time {1:F2} - Instrument type={2}, Number dependent scans={3}",
                            scan,
                            time,
                            scanDependents.RawFileInstrumentType,
                            scanDependents.ScanDependentDetailArray.Length);
                    }
                }
            }
        }
    }
}