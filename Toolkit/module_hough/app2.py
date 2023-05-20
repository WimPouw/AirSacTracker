import cv2
import numpy as np
import os

from tkinter import *


class App:
    def __init__(self):
        self.root = Tk()
        self.root.title("Video Processing")

        # Create the left window
        self.left_frame = Frame(self.root)
        self.left_frame.pack(side=LEFT)

        # Create the right window
        self.right_frame = Frame(self.root)
        self.right_frame.pack(side=RIGHT)

        # Create the video player widget
        self.video_player = cv2.VideoCapture()
        self.video_label = Label(self.left_frame)
        self.video_label.pack()

        # Create the process button
        self.process_button = Button(self.root, text="Process", command=self.process_video)
        self.process_button.pack()

        # Create the load button
        self.load_button = Button(self.root, text="Load", command=self.load_video)
        self.load_button.pack()

        # Create the output directory
        self.output_directory = "output"

        # Start the main loop
        self.root.mainloop()

    def load_video(self):
        # Get the video file name
        video_file = askopenfilename()

        # Load the video
        self.video_player.open(video_file)

        # Get the frame width and height
        width = int(self.video_player.get(3))
        height = int(self.video_player.get(4))

        # Create a blank frame for the original video
        original_frame = np.zeros((height, width, 3), np.uint8)

        # Show the original video in the left window
        while True:
            # Read the next frame
            ret, frame = self.video_player.read()

            # If the frame is not read, then the end of the video has been reached
            if not ret:
                break

            # Show the frame
            cv2.imshow("Left Window", frame)

            # Wait for a key press
            key = cv2.waitKey(0) & 0xFF

            # If the key pressed is `q`, then quit the program
            if key == ord("q"):
                cv2.destroyAllWindows()
                exit()

    def process_video(self):
        # Get the video file name
        video_file = askopenfilename()

        # Load the video
        self.video_player.open(video_file)

        # Get the frame width and height
        width = int(self.video_player.get(3))
        height = int(self.video_player.get(4))

        # Create an output video file
        output_file = os.path.join(self.output_directory, os.path.basename(video_file))
        out = cv2.VideoWriter(output_file, cv2.VideoWriter_fourcc(*'MJPG'), 25, (width, height))

        # Loop over the frames of the video
        while True:
            # Read the next frame
            ret, frame = self.video_player.read()

            # If the frame is not read, then the end of the video has been reached
            if not ret:
                break

            # Apply a simple processing function to the frame
            frame = cv2.blur(frame, (5, 5))

            # Write the processed frame to the output file
            out.write(frame)

        # Close the output file
        out.release()

        # Release the video capture object
        self.video_player.release()

        # Load the processed video file
        processed_video_file = os.path.join(self.output_directory, os.path.basename(video_file))
        processed_cap = cv2.VideoCapture(processed_video_file)

        # Get the frame width and height of the processed video
        processed_width = int(processed_cap.get(3))
        processed_height = int(processed_cap.get(4))

        # Create a blank frame for the processed video
        processed_frame = np.zeros((processed_height, processed_width, 3), np.uint8)

        # Show the processed video in the right window
        while True:
            # Read the next frame
            ret, frame = processed_cap.read()

            # If the frame is not read, then the end of the video has been reached
            if not ret:
                break

            # Show the frame
            cv2.imshow("Right Window", frame)

            # Wait for a key press
            key = cv2.waitKey(0) & 0xFF

            # If the key pressed is `q`, then quit the program
            if key == ord("q"):
                cv2.destroyAllWindows()
                exit()


if __name__ == "__main__":
    app = App()