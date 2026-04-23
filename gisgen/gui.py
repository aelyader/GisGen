import os

from PySide6.QtCore import QTime, QTimer, Qt
from PySide6.QtGui import QFont, QPixmap
from PySide6.QtWidgets import (
    QFileDialog,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QMessageBox,
    QProgressBar,
    QPushButton,
    QStyle,
    QTextEdit,
    QVBoxLayout,
    QWidget,
    QCheckBox,
)

from .worker import WorkerThread


class ProgressWindow(QWidget):
    def __init__(self):
        super().__init__()
        # Declare instance attributes here
        self.startTime = None
        self.timer = QTimer(self)
        self.layout = QVBoxLayout(self)
        self.elapsedTimeLabel = None  # Declare now, initialize in initUI
        self.main_progress_bar = None  # Declare now, initialize in initUI
        self.text_box = None  # Declare now, initialize in initUI
        # Connect the timer's timeout signal to updateElapsedTime
        self.timer.timeout.connect(self.updateElapsedTime)
        self.initUI()

    def initUI(self):
        self.setWindowTitle("Progress")

        # Elapsed time label
        self.elapsedTimeLabel = QLabel("00:00", self)
        self.layout.addWidget(self.elapsedTimeLabel)

        self.main_progress_bar = QProgressBar(self)
        self.layout.addWidget(self.main_progress_bar)

        self.text_box = QTextEdit(self)
        self.text_box.setReadOnly(True)
        self.layout.addWidget(self.text_box)

        self.setGeometry(500, 500, 600, 600)

    def startElapsedTimeTimer(self):
        self.startTime = QTime.currentTime()  # Record the start time
        self.timer.start(1000)  # Start the timer to update every second

    def updateElapsedTime(self):
        if self.startTime:
            # Calculate the elapsed time
            now = QTime.currentTime()
            elapsed_time = self.startTime.secsTo(now)

            # Convert seconds to minutes and seconds
            minutes, seconds = divmod(elapsed_time, 60)
            # Format the string to display
            time_string = "{:02d}:{:02d}".format(minutes, seconds)

            self.elapsedTimeLabel.setText(time_string)

    def update_progress(self, value):
        self.main_progress_bar.setValue(value)

    def update_message(self, message):
        # Check if the message starts with any of the specified keywords
        if message.startswith("Step"):
            # Apply HTML formatting for bold and larger font, starting with a line break
            formatted_message = f"<br><div style='font-weight:bold; font-size:14pt; margin-left: 0px;'>{message}</div>"
        elif message.startswith(("->", "Skipping", "Mismatch")):
            # For lines starting with "->", add indentation using a div for proper margin application
            formatted_message = f"<br><div style='font-size:12pt; margin-left: 35px;'>{message}</div>"
        elif message.startswith(("MD5", "Downstream")):
            formatted_message = f"<br><div style='margin-left: 0px; font-weight:bold; color: green; font-size:18pt; margin-left: 50px;'>{message}</div><br>"
        else:
            # Apply a smaller font size increase for other messages, starting with a line break
            formatted_message = f"<br><div style='font-size:12pt; margin-left: 0px;'>{message}</div>"

        # Append the formatted message to the text box
        self.text_box.append(formatted_message)


class MainApplication(QWidget):
    output_folder = ""

    def __init__(self):
        super().__init__()
        self.worker_thread = WorkerThread()
        self.main_layout = QVBoxLayout()
        self.perform_downstream_analysis = False
        self.worker_thread.resume_signal.connect(self.worker_thread.perform_heavy_tasks)
        self.gisaid_entry = QLineEdit()
        self.genbank_entry = QLineEdit()  # Assuming similar pattern for other QLineEdit attributes
        self.gisaid_fasta_folder_entry = QLineEdit()
        self.output_folder_entry = QLineEdit()
        self.filter_omicron_checkbox = QCheckBox()
        self.downstream_analysis_checkbox = QCheckBox()
        self.progress_window = None
        self.initUI()

    def initUI(self):
        self.setWindowTitle("GISAID <--> GenBanks Sample Matching")
        self.setLayout(self.main_layout)
        self.addComponents()
        self.worker_thread.pause_execution.connect(self.wait_for_user_input)

    def addComponents(self):
        self.main_layout.addLayout(self.createSubtitleLayout())
        self.main_layout.addLayout(self.createGISAIDLayout())
        self.main_layout.addLayout(self.createGenBankLayout())
        self.main_layout.addLayout(self.createOutputFolderLayout())
        self.main_layout.addLayout(self.createGISAIDFASTAFolderLayout())
        self.main_layout.addLayout(self.createFilterOmicronLayout())
        self.main_layout.addLayout(self.createDownstreamAnalysisCheckbox())
        self.main_layout.addLayout(self.createButtonsLayout())

    def createSubtitleLayout(self):
        # Create a vertical layout
        vertical_layout = QVBoxLayout()

        # Create a label for the logo and set the pixmap
        logo_label = QLabel(self)
        # Get the directory in which the script is located
        package_dir = os.path.dirname(os.path.abspath(__file__))
        project_root = os.path.dirname(package_dir)
        image_path = os.path.join(project_root, "images", "Final_gisgen_logo_new.jpg")
        logo_pixmap = QPixmap(image_path)  # Adjust path if necessary
        if logo_pixmap.isNull():  # Check if the QPixmap failed to load
            print("The logo image was not found. Displaying placeholder text instead.")
            logo_label.setText("Logo Not Found")  # Placeholder text
        else:
            logo_label.setPixmap(logo_pixmap.scaled(250, 250, Qt.KeepAspectRatio))

        vertical_layout.addWidget(logo_label, alignment=Qt.AlignCenter)

        # Create and configure the Info Help button
        info_button = QPushButton("Info Help \n(Step by Step Guide)\n")
        font = QFont()
        font.setPointSize(11)
        info_button.setFont(font)
        info_button.clicked.connect(self.show_instructions)
        vertical_layout.addWidget(info_button)  # Add the button below the logo

        return vertical_layout

    def createGISAIDLayout(self):
        gisaid_layout = QHBoxLayout()
        gisaid_label_text = """<b><span style="font-size:12pt;">GISAID Meta Table:</span></b><br>
            [ .tar.xz sample metatable from <br><a href="https://gisaid.org/">https://gisaid.org/</a> ]"""
        gisaid_label = QLabel(gisaid_label_text)

        gisaid_button = QPushButton("Browse")
        gisaid_button.clicked.connect(lambda: self.select_file(self.gisaid_entry))
        gisaid_layout.addWidget(gisaid_label)
        gisaid_layout.addWidget(self.gisaid_entry)
        gisaid_layout.addWidget(gisaid_button)
        return gisaid_layout

    def createGenBankLayout(self):
        genbank_layout = QHBoxLayout()
        genbank_label = QLabel(
            """<b><span style="font-size:12pt;">GenBank Meta Table:</span></b><br> [ CSV sample metatable from <br> <a href="https://www.ncbi.nlm.nih.gov/labs/virus">https://www.ncbi.nlm.nih.gov/labs/virus</a> ]""")
        genbank_button = QPushButton("Browse")
        genbank_button.clicked.connect(lambda: self.select_file(self.genbank_entry))
        genbank_layout.addWidget(genbank_label)
        genbank_layout.addWidget(self.genbank_entry)
        genbank_layout.addWidget(genbank_button)
        return genbank_layout

    def createOutputFolderLayout(self):
        output_folder_layout = QHBoxLayout()
        output_folder_label = QLabel(
            """<b><span style="font-size:12pt;">Output Folder:</span></b><br> [ Select a folder where results will be saved ]""")
        output_folder_button = QPushButton("Browse")
        output_folder_button.clicked.connect(lambda: self.select_folder(self.output_folder_entry))
        output_folder_layout.addWidget(output_folder_label)
        output_folder_layout.addWidget(self.output_folder_entry)
        output_folder_layout.addWidget(output_folder_button)
        return output_folder_layout

    def createGISAIDFASTAFolderLayout(self):
        gisaid_fasta_folder_layout = QHBoxLayout()
        gisaid_fasta_folder_label = QLabel(
            """<b><span style="font-size:12pt;">GISAID FASTA Folder:</span></b><br> [ Select a folder for the GISAID FASTA.<br> <span style="color:red;">This folder needs to be empty.</span> ]<br>""")
        self.gisaid_fasta_folder_entry = QLineEdit()
        gisaid_fasta_folder_button = QPushButton("Browse")
        gisaid_fasta_folder_button.clicked.connect(lambda: self.select_folder(self.gisaid_fasta_folder_entry))
        gisaid_fasta_folder_layout.addWidget(gisaid_fasta_folder_label)
        gisaid_fasta_folder_layout.addWidget(self.gisaid_fasta_folder_entry)
        gisaid_fasta_folder_layout.addWidget(gisaid_fasta_folder_button)
        return gisaid_fasta_folder_layout

    def createFilterOmicronLayout(self):
        filter_omicron_layout = QHBoxLayout()
        filter_omicron_layout.addStretch(1)
        self.filter_omicron_checkbox = QCheckBox("[ Filter by Omicron ]")  # Define the checkbox as an attribute
        font = QFont()
        font.setPointSize(11)  # Adjust the size as needed
        self.filter_omicron_checkbox.setFont(font)

        filter_omicron_layout.addWidget(self.filter_omicron_checkbox)
        filter_omicron_layout.addStretch(1)
        return filter_omicron_layout

    def createDownstreamAnalysisCheckbox(self):
        downstream_analysis_layout = QHBoxLayout()
        downstream_analysis_layout.addStretch(1)

        self.downstream_analysis_checkbox = QCheckBox("[ Enable Downstream Analysis ]")
        font = QFont()
        font.setPointSize(11)  # Adjust the size as needed
        self.downstream_analysis_checkbox.setFont(font)

        downstream_analysis_layout.addWidget(self.downstream_analysis_checkbox)
        analysis_info_button = QPushButton()
        analysis_info_button.setIcon(self.style().standardIcon(getattr(QStyle, "SP_MessageBoxInformation")))
        analysis_info_button.clicked.connect(self.show_analysis_description)  # Connect to the new method

        downstream_analysis_layout.addWidget(analysis_info_button)
        downstream_analysis_layout.addStretch(1)
        return downstream_analysis_layout

    def createButtonsLayout(self):
        buttons_layout = QHBoxLayout()
        run_button = QPushButton("Run")
        font = QFont()
        font.setPointSize(14)  # Adjust the size as needed
        font.setBold(True)  # Make the text bold
        run_button.setFont(font)

        run_button.clicked.connect(self.run_script)
        buttons_layout.addWidget(run_button)
        return buttons_layout

    def run_script(self):
        # Warning message about computational power and GUI suspension
        warning_msg = QMessageBox()
        warning_msg.setIcon(QMessageBox.Icon.Warning)
        warning_msg.setWindowTitle("Warning")
        warning_msg.setText(
            "<p><span style='font-weight:bold; font-size:13pt;'>WARNING</span></p>"
            "<p><span font-size:13pt;'>This program needs substantial computational power and may temporarily suspend the GUI to allocate more threads for processing, depending on your hardware's capabilities.</span></p>"
            "<p><span style='font-weight:bold; font-size:13pt;'>**** If the GUI seems unresponsive, please avoid closing the application. ****</span></p>")
        warning_msg.setStandardButtons(
            QMessageBox.StandardButton.Ok | QMessageBox.StandardButton.Cancel)  # Uncomment for PyQt6

        response = warning_msg.exec()

        if response == QMessageBox.StandardButton.Ok:
            if not self.worker_thread.isRunning():
                # Set up the variables for the worker thread
                self.worker_thread.gisaid_data_file = self.gisaid_entry.text()
                self.worker_thread.genbank_data_file = self.genbank_entry.text()
                self.worker_thread.filter_by_omicron = self.filter_omicron_checkbox.isChecked()
                self.worker_thread.output_folder = self.output_folder_entry.text()
                self.worker_thread.gisaid_fasta_folder_path = self.gisaid_fasta_folder_entry.text()
                self.worker_thread.perform_downstream_analysis = self.downstream_analysis_checkbox.isChecked()
                MainApplication.output_folder = self.output_folder_entry.text()
                self.worker_thread.output_folder = MainApplication.output_folder

                # Start the worker thread
                self.worker_thread.start()

                # Initialize and display the progress window
                self.progress_window = ProgressWindow()  # Instantiate ProgressWindow when needed
                self.progress_window.show()
                self.worker_thread.progress_update.connect(self.progress_window.update_progress)
                self.worker_thread.update_message.connect(self.progress_window.update_message)
                self.progress_window.startElapsedTimeTimer()

        elif response == QMessageBox.StandardButton.Cancel:
            # If the user cancels, do not start the process
            pass

    def resume_script(self):
        if self.downstream_analysis_checkbox.isChecked():
            self.worker_thread.perform_downstream_analysis = True
        self.worker_thread.resume_work()

    def select_file(self, line_edit):
        file_path, _ = QFileDialog.getOpenFileName(self, "Select File")
        if file_path:  # This checks if a file was actually selected
            line_edit.setText(file_path)

    def select_folder(self, line_edit):
        folder_path = QFileDialog.getExistingDirectory(self, "Select Folder")
        if folder_path:  # This checks if a folder was actually selected
            line_edit.setText(folder_path)

    def wait_for_user_input(self):
        # Create a modal dialog box to wait for user input
        gisaid_fasta_path = self.worker_thread.gisaid_fasta_folder_path
        file_path = self.worker_thread.sample_ids_file_path
        msgBox = QMessageBox()
        msgBox.setWindowTitle("Pause Execution")
        msgBox.setText(f"\n         **********   STOP   **********\n"
                       f"\nThe Code has PAUSED. Please complete the required steps and then press 'Continue'."
                       f"\nPlease Follow the steps below:\n"
                       f"\nStep 1:"
                       f"\n ----> Open a web browser and Go to https://www.epicov.org/ "
                       f"\n     ----> login using your GISAID credentials\n"
                       f"\nStep 2:"
                       f"\n ----> In https://www.epicov.org/ select [EpiCoV] and under EpiCoV, select [Search]"
                       f"\n     ----> At the bottom right hand side between EPI_SET and Analysis box, click [Select]"
                       f"\n         ----> Click [Choose File] and select the following text file:"
                       f"\n             ---->  {file_path}\n"
                       f"\nStep 3:"
                       f"\n ----> Download multi-fasta file for the queried sample set\n"
                       f"\nStep 4:"
                       f"\n ----> Place the single multi-fasta file in the following folder:  "
                       f"\n     ----> {gisaid_fasta_path}"
                       f"\n         ----> Make sure the single multi-fasta file is the only file in the folder\n"
                       f"\n Now you can Click 'OK' to proceed.\n")

        msgBox.setStandardButtons(QMessageBox.StandardButton.Ok)
        if msgBox.exec() == QMessageBox.StandardButton.Ok:
            # Logic after the user presses 'Continue' on the dialog box
            # This will resume the worker thread
            self.worker_thread.resume_work()

    def show_analysis_description(self):
        description = """
        Downstream Analysis Description:
        --------------------------------
        This analysis includes:
        - Temporal Analysis of Variants: 
            Tracks the evolution of virus variants over time.

        - Geographical Distribution: 
            Studies the spread of variants across different regions.

        - Sequence Length Variation: 
            Analyzes the variation in sequence lengths among different samples.

        - Amino Acid Substitutions: 
            Identifies common substitutions in the virus sequences.

        - Submission Date Lag: 
            Examines the time lag between sample collection and submission dates.

        Enable this option to include these analyses in the final report.
        """
        QMessageBox.information(self, "Downstream Analysis Description", description)

    def show_instructions(self):
        instructions = """
        1. Provide GISAID Meta Table:  
            Click the "Browse" button next to "GISAID Meta Table" and 
            select the .tar.xz file you downloaded from https://gisaid.org/.

        2. Provide GenBank Meta Table:  
            Click the "Browse" button next to "GenBank Meta Table" and 
            select the CSV file you downloaded from 
            https://www.ncbi.nlm.nih.gov/labs/virus.

        3. Select Output Folder:  
            Select a folder where the matching table and 
            other result output files will be saved.

        4. GISAID FASTA Folder:  
            Choose an empty folder for the GISAID multi-sample FASTA file. 
            Ensure it's the only file in the folder when prompted to add it.

        5. Filter Options:  
            Optionally, you can select the checkboxes to filter by Omicron variant.

        5. Enable Downstream Analysis:  
            You can select the checkboxes to perform several downstream analysis: 
            [Temporal Analysis of Variants] [Geographical Distribution] 
            [Sequence Length Variation] [Amino Acid Substitutions] 
            [Submission Date Lag]

        6. Run:  
            Click the "Run" button to execute the script.

        7. Monitor Progress:  
            A second window with a progress bar and task description will appear. 
            This program demands substantial computational power and may briefly 
            suspend the GUI to dedicate more threads to processing, depending on 
            your hardware. 

        ** If the GUI becomes unresponsive, please do not close the application **

        """
        QMessageBox.information(self, "Instructions", instructions)
