	var topNavNew = ['Home','News','Missions','Multimedia','Connect','Aboutus'];
function addMenuSupport() {
    // create, and post process menus

    var headerstring = '<div style="background-color: #1f1f1f; height: 100px;">' +
        '<h3><a href="http://www.nasa.gov"><div style="float: left;">' +
        '<img src="http://www.nasa.gov/externalflash/cmsdev_css/header_logo_revised.gif" alt="National Aeronautics and Space Administration">' +
        '</a></div><div style="color:#FFF;"><br>' +
        '<b>Disclaimer</b>: This material is being kept online for historical purposes. ' +
        'Though accurate at the time of publication, it is no longer being updated. ' +
        'The page may contain broken links or outdated information, ' + 
        'and parts may not function in current web browsers. ' + 
        'Visit <a href="http://www.nasa.gov">NASA.gov</a> for current information';

    document.getElementById('top_header').innerHTML = headerstring;
}

		function createMajorNav(menuParent, menuSrc, menuClass) {
				       // set up some vars
					//var logoItem  = logoItem;
					var navParent = menuParent;
					var navData	=	menuSrc;
					// create menu root dom
					/*var navBucket            = document.createElement('div');
						navBucket.id         = menuId;
						navBucket.className	 = menuClass;*/
					var navBucketDiv		 =	document.createElement('div');
						navBucketDiv.id		 =	menuClass;
				// loop through level 1
					
					//for (each in navData) {
						navData.each(function(value,index){
						var topLevel = navData[index];
						var topName  = topLevel[0];
						var topLink  = topLevel[1];
                                       
                                  
                                   
					              // create first level item DOM
                                                 var navItemA;
                                                 var navItemASpan;
                                                 var navItem  = document.createElement('h2');
                                                 if(topName.toLowerCase()=='connect') {
                                                 navItem.style.width="125px";
                                                 }
							navItem.id		   = 'nav-'+topName.replace(" ","").toLowerCase();
							if(topName.toLowerCase()=='home'){
								navItem.className	=	'nav-'+topName.replace(" ","").toLowerCase();
							}
                                                 if(topName.toLowerCase()=='connect') 
                                                 {
                                                 navItemA           = document.createElement('a');
						       navItemA.href      = topLink;
                                                 navItemA.style.width="125px";
                                                 navItemA.style.backgroundImage="url(/templateimages/redesign/navigation/TopNav/navtitle-connect.gif)";
                                                                                                  
                                                 navItemASpan			=	document.createElement('span');
                                                 navItemASpan.style.display="none";
							navItemASpan.innerHTML = topName;
                                                 }
                                                 else
                                                 {
                                                 navItemA = document.createElement('a');
							navItemA.href = topLink;
							navItemASpan =document.createElement('span');
							navItemASpan.innerHTML = topName;
                                                 }
						
						navItemA.appendChild(navItemASpan);
						navItem.appendChild(navItemA);
						
				
						if (topLevel.length > 2) {
							var subLevel  = topLevel[2];
						// loop through level 2
							//for (each in subLevel) {
								var subBucket = document.createElement('div');
								if(topName.toLowerCase()=='home'){
									subBucket.className	=	"dropper 129";
								}
                                                        else if (topName.toLowerCase()=='connect'){
									subBucket.className	=	"dropper 125";
								}
                                                        else{
									subBucket.className	=	"dropper 140";
								}
								var subBucketUl			=	document.createElement('ul');

                                                  
                                                        if(topName.toLowerCase()=='connect'){
									subLevel.each(function(value,index){

								var subName = subLevel[index][0];
								var subLink = subLevel[index][1];

                                               
							 //alert("subname " +	subName+ "sublink "+subLink);
                                                  // create sublevel item element
								var subItem            = document.createElement('li');
                                                              subItem.style.lineHeight="18px";
                                                              subItem.style.height="20px";

							 
                                                        //alert(subLevel[index][2]);
								if (subLevel[index][2] && subLevel[index][2] != null && subLevel[index][2] != eval(""))
                                                           {
                                                            var subImage = subLevel[index][2];
                                                            var subItemB = document.createElement('a');
                                                        	subItemB.href      = subLink;
									subItemB.innerHTML = subImage;
                                                               subItem.appendChild(subItemB);
                                                           }
	

                                                               var subItemA           = document.createElement('a');
                                                        	subItemA.href      = subLink;
									subItemA.innerHTML = subName;
								
                                                           
                                                        subItem.appendChild(subItemA);
                                                 	subBucketUl.appendChild(subItem);
								
						
							});
							}
                       


								if(topName.toLowerCase()!='home' && topName.toLowerCase()!='connect'){
									subLevel.each(function(value,index){

								var subName = subLevel[index][0];
								var subLink = subLevel[index][1];
							//alert("subname " +	subName+ "sublink "+subLink);
							// create sublevel item element
								var subItem            = document.createElement('li');
								var subItemA           = document.createElement('a');
									subItemA.href      = subLink;
									subItemA.innerHTML = subName;
									
								subItem.appendChild(subItemA);
								subBucketUl.appendChild(subItem);
								
						
							});
							}
							subBucket.appendChild(subBucketUl);
							
						}
						
							navBucketDiv.appendChild(navItem);					
							navBucketDiv.appendChild(subBucket);
						
                                   
					});

					navParent.appendChild(navBucketDiv);
					
				}
				
				
	function addTopNavContent() {

	var navData = topNavNew;

	var bodyTag= document.getElementsByTagName('body')[0];

	var mainNavTag = document.getElementsByTagName('main-nav');

	var i;
        
	var mainDivTag = document.getElementById('main-nav');
	//														alert('maintag :' + mainDivTag);
		var subDivTag = document.createElement('div');
		subDivTag.id = "dropper_wrapper";

	for(i=0;i<navData.length;i++)
	{
	
		var headerTag = document.createElement('h2');
		if(navData[i].toLowerCase() == 'aboutus')
		{
			navData[i] = "aboutnasa";
		}
		headerTag.id="shelf-nav-"+navData[i].toLowerCase();
		headerTag.className="nav-"+navData[i].toLowerCase()+"-out";
		
		if(navData[i].toLowerCase() == 'home') {
				headerTag.style.width="129px"
			}else if(navData[i].toLowerCase() == 'connect') {
				headerTag.style.width="124px"
			}
			var anchorTag = document.createElement('a');
			if(navData[i].toLowerCase() == 'aboutnasa')
			{
				navData[i] = "aboutus";
			}
			anchorTag.id = "topnav_"+navData[i].toLowerCase()+"link";
			anchorTag.rel="topnav_"+navData[i].toLowerCase()+"link_submenu";
			if(navData[i].toLowerCase() == 'home') {
				anchorTag.style.width="129px"
			}else if(navData[i].toLowerCase() == 'connect') {
				anchorTag.style.width="124px"
			}
			if(navData[i].toLowerCase() == 'home') {
				anchorTag.href="/home/index.html";
			}else if(navData[i].toLowerCase() == 'news') {
				anchorTag.href="/news/index.html";
			}else if(navData[i].toLowerCase() == 'missions') {
				anchorTag.href="/missions/index.html";
			}else if(navData[i].toLowerCase() == 'multimedia') {
				anchorTag.href="/multimedia/index.html";
			}else if(navData[i].toLowerCase() == 'connect') {
				anchorTag.href="/connect/index.html";
			}else {
				anchorTag.href="/about/index.html";
			}

			if(navData[i].toLowerCase() != 'home')
			{
				anchorTag.rev="/templateimages/redesign/shelfnav/"+navData[i].toLowerCase()+"_topnav.html";
			}
				var spanTag = document.createElement('span');
				if(navData[i].toLowerCase() == 'connect')
				{
					//spanTag.setAttribute("style","display:none;");
					//spanTag.style.display="none";
				}
				if((navData[i].replace(" ","")).toLowerCase() == 'aboutus')
				{
					navData[i] = "About Us";
				}
				spanTag.style.display="none";
			    spanTag.innerHTML = navData[i];
			anchorTag.appendChild(spanTag);
			headerTag.appendChild(anchorTag);
			
			var divTag = document.createElement("div");
			divTag.id = "topnav_"+(navData[i].replace(" ","")).toLowerCase()+"link_submenu";			 
			divTag.className ="dropdownmenu";			 
			
			
			subDivTag.appendChild(headerTag);
			subDivTag.appendChild(divTag);	
			}

			mainDivTag.appendChild(subDivTag);

		//	alert('done :');

		//dropdowncontent.init("topnav_homelink", "right-bottom", 100, "mouseover")				  
		dropdowncontent.init("topnav_newslink", "right-bottom", 700, "mouseover")
		dropdowncontent.init("topnav_missionslink", "right-bottom", 700, "mouseover")
		dropdowncontent.init("topnav_multimedialink", "right-bottom", 700, "mouseover")
		dropdowncontent.init("topnav_connectlink", "right-bottom", 700, "mouseover")
		dropdowncontent.init("topnav_aboutuslink", "right-bottom", 700, "mouseover")

	}
			




			
