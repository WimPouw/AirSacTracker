import cv2
import tkinter as tk
from tkinter import filedialog

class VideoProcessorGUI:

    def __init__(self):
        self.window = tk.Tk()
        self.window.title("Video Processor")

        # Create buttons
        self.load_button = tk.Button(self.window, text="Load Video", command=self.load_video)
        self.process_button = tk.Button(self.window, text="Process Video", command=self.process_video)
        self.quit_button = tk.Button(self.window, text="Quit", command=self.window.quit)

        # Create video panels
        self.original_video_panel = tk.Label(self.window)
        self.processed_video_panel = tk.Label(self.window)

        # Pack buttons and video panels
        self.load_button.pack(side="top", padx=10, pady=10)
        self.process_button.pack(side="top", padx=10, pady=10)
        self.quit_button.pack(side="top", padx=10, pady=10)
        self.original_video_panel.pack(side="left", padx=10, pady=10)
        self.processed_video_panel.pack(side="right", padx=10, pady=10)

        # Initialize video variables
        self.original_video = None
        self.processed_video = None

    def load_video(self):
        # Open file dialog to select video file
        file_path = filedialog.askopenfilename(filetypes=[("Video Files", "*.mp4")])

        # Load video and display it on the original video panel
        self.original_video = cv2.VideoCapture(file_path)
        ret, frame = self.original_video.read()
        frame = cv2.cvtColor(frame, cv2.COLOR_BGR2RGB)
        img = tk.PhotoImage(width=frame.shape[1], height=frame.shape[0])
        img.blank()
        img.paste(tk.PhotoImage(image=tk.Image.fromarray(frame)))
        self.original_video_panel.configure(image=img)
        self.original_video_panel.image = img

    def process_video(self):
        if self.original_video is None:
            return

        # Process video and save it to file
        output_file_path = "processed_video.mp4"
        output_video = cv2.VideoWriter(output_file_path, cv2.VideoWriter_fourcc(*'mp4v'), 30, (640, 480))
        while True:
            ret, frame = self.original_video.read()
            if not ret:
                break
            # Apply some processing to the frame
            processed_frame = cv2.flip(frame, 1)
            # Write the processed frame to the output video
            output_video.write(processed_frame)

        output_video.release()

        # Load the processed video and display it on the processed video panel
        self.processed_video = cv2.VideoCapture(output_file_path)
        ret, frame = self.processed_video.read()
        frame = cv2.cvtColor(frame, cv2.COLOR_BGR2RGB)
        img = tk.PhotoImage(width=frame.shape[1], height=frame.shape[0])
        img.blank()
        img.paste(tk.PhotoImage(image=tk.Image.fromarray(frame)))
        self.processed_video_panel.configure(image=img)
        self.processed_video_panel.image = img

    def run(self):
        self.window.mainloop()

if __name__ == '__main__':
    gui = VideoProcessorGUI()
    gui.run()