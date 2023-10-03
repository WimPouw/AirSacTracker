<h1 align="center" style=font-size:200px>A toolkit for the dynamic study of elastic kinematics.</h1>

<video src='https://tsg-131-174-75-200.hosting.ru.nl/samples_airsactoolkit/June16_02_circle_rec.mp4' width=180/> 

Biological structures are defined by elements like bones and cartilage, and elastic elements like muscles and membranes. Computer vision advances have enabled automatic tracking of animal skeletal poses. However, the elastic and soft-tissues of organisms, like the nose of Elephant seals, or the buccal sac of frogs, have been poorly studied as no computer vision methods are optimized for tracking such elastic kinematics. This leaves major gaps in different areas in biology. In the area of primatology, most critically, the function of air sacs is widely debated and many questions exist about their role in communication and human language evolution. Moving towards the dynamic study of soft-tissue elastic structures, we present a toolkit for the automated tracking of semi-circular elastic structures in biological video data. The toolkit contains unsupervised computer vision tools (using Hough transform) and supervised deep learning (by adapting Deeplabcut) methodology to track inflation of laryngeal air sacs or other biological spherical objects (e.g., gular cavities). 

# Toolkit

<table>
  <thead>
    <tr>
      <th>✨</th>
      <th>Feature</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>✅</td>
      <td>Hough transform to detect semi-circles (unsupervised method): 
https://wimpouw.github.io/AirSacTracker/AirSacTracking_with_Hough.html</td>
    </tr>
 <tbody>
    <tr>
      <td>✅</td>
      <td>Deeplabcut + Landau circle estimation (supervised method): 
https://wimpouw.github.io/AirSacTracker/DLC+.html</td>
    </tr>
      <td>✅</td>
      <td>Open dataset on Siamang Air Sacs: 
https://data.donders.ru.nl/collections/mine?3</td>
    </tr>
</table>

## Installation / requirements ## 
See requirements.txt for each module.

### Pipeline ### 

<img src = ...  > (pipeline image) 

### file structure ###

--Project -> contains all the code and materials for the manuscript
--Docs -> contains all the github pages
 
--Toolkit/
	--Input/
		--inputvideos.mp4
	--Module_DLC+/
		--DLC+.ipynb
		--DLC/
			--Meta information/
		--Output_DLC/
		  -- DLC_labeled_videos.mp4 
		  -- DLC_coordinates.csv
		  -- DLC_coordinates.h5
		  -- inputvideos.mp4
		--Output_DLC+/
			--timeseries/
				--circlecoordinates.csv
			--DLC+_labeled_videos.mp4
	--Module_Hough/
		-- AirSacTracking_with_Hough.ipynb
		-- results/
			-- hough_labeled_videos.mp4 
			-- coordinates.csv
			

# Code contributers
Lara Burchardt, Yana van der Sande, Wim Pouw

##Reference ## 
TBA




