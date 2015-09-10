awk '{newlon=$2; if(newlon< 130) newlon=360+newlon;}{print $1, newlon, $3}' haz_pts_W-65.txt  > tmp.txt
head -n1 haz_pts_W-65.txt > haz_pts_W130.txt
awk 'NR>1' tmp.txt >> haz_pts_W130.txt
