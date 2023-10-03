<h1 align="center" style=font-size:200px>A toolkit for the dynamic study of elastic kinematics</h1>
<h2 align="center" style=font-size:200px>Burchardt (burchardt@leibniz-zas.de), van der Sande, Kehy, Gamba, Ravignani, Pouw (wim.pouw@donders.ru.nl)</h2>

<a name="overview"></a>

<img src = /docs/videos/side_by_side.gif  >

Biological structures are defined by elements like bones and cartilage, and elastic elements like muscles and membranes. Computer vision advances have enabled automatic tracking of animal skeletal poses. However, the elastic and soft-tissues of organisms, like the nose of Elephant seals, or the buccal sac of frogs, have been poorly studied as no computer vision methods are optimized for tracking such elastic kinematics. This leaves major gaps in different areas in biology. In the area of primatology, most critically, the function of air sacs is widely debated and many questions exist about their role in communication and human language evolution. Moving towards the dynamic study of soft-tissue elastic structures, we present a toolkit for the automated tracking of semi-circular elastic structures in biological video data. The toolkit contains unsupervised computer vision tools (using Hough transform) and supervised deep learning (by adapting Deeplabcut) methodology to track inflation of laryngeal air sacs or other biological spherical objects (e.g., gular cavities). 

# Toolkit

<table>
  <thead>
    <tr>
      <th></th>
      <th>Component</th>
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

## Manuscript pipeline

<img src = /docs/images/Workflow_figure.png  > (pipeline image) 


## Installation / requirements ## 
See requirements.txt for each module. You can install the requirements by entering in your terminal 'pip -r requirements.txt' (after navigating to the folder where the requirements.txt is located)

## File structure

- Project -> contains all the code and materials for the manuscript
- Docs -> contains all the github pages
- Toolkit -> contains the hough and DLC+ code
	
## Code contributers
Lara Burchardt, Yana van der Sande, Wim Pouw

## Reference
TBA




