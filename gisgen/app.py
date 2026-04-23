import cProfile
import os
import sys

from PySide6.QtWidgets import QApplication

from .gui import MainApplication


def main():
    app = QApplication(sys.argv)
    window = MainApplication()
    window.show()

    with cProfile.Profile() as profiler:
        app.exec()

    profile_output_path = os.path.join(MainApplication.output_folder, "profile.prof")
    profiler.dump_stats(profile_output_path)


if __name__ == "__main__":
    main()
