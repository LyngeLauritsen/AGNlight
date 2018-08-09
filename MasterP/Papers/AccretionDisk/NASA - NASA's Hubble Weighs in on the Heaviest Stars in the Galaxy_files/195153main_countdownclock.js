/**
* Utility function for displaying a up/down time counter.
*
* @param inputDate :: The date of the launch/mission
* @param imgSuffix :: The suffix used for the images
* @param imgPath :: The relative path for the images
*/
function MissionTimer(inputDate, divId) {

	//added arrays to replace the time zones for opera browser issue
	var inputString = inputDate.toString();
	var replacedString = null;
	var timeZone =        new Array("EDT",      "EST",      "PDT",      "PST" ,      "CDT",      "CST",      "MDT",      "MST",      "AKDT",     "AKST",     "ADT",      "AST",      "HST");
	var timeZoneReplace = new Array("UTC-0400", "UTC-0500", "UTC-0700", "UTC-0800" , "UTC-0500", "UTC-0600", "UTC-0600", "UTC-0700", "UTC-0800", "UTC-0900", "UTC-0300", "UTC-0400", "UTC-1000");

	for(var i=0;i<timeZone.length;i++){
		if(inputString.match(timeZone[i])==timeZone[i]){
		  replacedString = inputString.replace(timeZone[i],timeZoneReplace[i]);
		}
	}

	var now = new Date();
	gmt_now = now.toUTCString();   
	now_ms = Date.parse(gmt_now)
	var launch = new Date(replacedString);
	var gmt_launch = launch;
	gmt_launch = launch.toUTCString();
	launch_ms = Date.parse(gmt_launch);
	var gap =launch_ms-now_ms;

	if (gap < 0) {
		gap = Math.abs(gap);
	}


	var day_gap_raw = (gap/(1000*60*60*24));
	var hr_gap_raw = (gap/(1000*60*60));
	var min_gap_raw = (gap/(1000*60));

	var day_gap = Math.floor(gap/(1000*60*60*24));
	var hr_gap = Math.floor(hr_gap_raw-(day_gap * 24));     
	// changed this to Math.floor from Math.round to see if it fixes the "final minute" problem
	var mn_gap = Math.floor((min_gap_raw-(day_gap*24*60))-(hr_gap*60));

	// Calculate the number of seconds left after minutes are calculated.
	var min_gap_floor = Math.floor(min_gap_raw);
	sec_gap = Math.round((min_gap_raw - min_gap_floor) *60);

	var daytho = Math.floor(day_gap/1000);
	if (daytho < 1 ) {
		daytho = 0;
	}

	var dayhun = Math.floor(day_gap / 100);
	var dayten = Math.floor((day_gap - (dayhun * 100))/10);
	var dayten_raw = (day_gap - (dayhun * 100))/10;
	var dayone = Math.floor((dayten_raw*10) - (dayten*10));

	var dayhunabs = Math.abs(dayhun);

	if (dayhunabs >= 10) {
		dayhun_div = Math.floor(dayhun/10);
		dayhun = dayhun - (dayhun_div * 10);
	}

	var hrten = Math.floor((hr_gap)/10);
	var hrone = Math.floor((hr_gap) - (hrten *10));

	var mnten = Math.floor((mn_gap)/10);
	var mnone = Math.floor((mn_gap) - (mnten * 10));

	var secten = Math.floor(sec_gap/10);
	var secone = Math.floor(sec_gap - (secten*10));
	
	var day = dayone;
	if (dayten > 0) day = ''+dayten+dayone;
	if (dayhun > 0) day = ''+dayhun+dayten+dayone;
	if (daytho > 0) day = ''+daytho+dayhun+dayten+dayone;
	

        var htmlSnippet = '<div id="day">'+day
        		+ '</div><div id="hour">'+hrten+hrone 
        		+ '</div><div id="minute">'+mnten+mnone
        		+ '</div><div id="second">'+secten+secone
        		+ '</div>';
        		
        		
	document.getElementById(divId).innerHTML = htmlSnippet;

	// recursive call to the function on every second
	setTimeout("MissionTimer('" + inputDate + "', '" + divId + "')", 1000);
}







