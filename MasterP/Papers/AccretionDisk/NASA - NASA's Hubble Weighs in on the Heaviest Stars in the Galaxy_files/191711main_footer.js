function createFooterlogo(){
				var footerParent	=	document.getElementById('footer');
				var logo = createLogo('NASA Home', 'nasa_logo_footer');
				var footerp			=	 document.createElement('p');
					footerp.id		=	"footercol1";
					
					var footer2	=	createFooter(document.getElementById('footer'),footercol2,'footercol2');
					var footer3	=	createFooter(document.getElementById('footer'),footercol3,'footercol3');
					var footer4	=	createFooter(document.getElementById('footer'),footercol4,'footercol4');

					footerParent.appendChild(logo);
					footerParent.appendChild(footerp);
					footerParent.appendChild(footer2);
					footerParent.appendChild(footer3);
					footerParent.appendChild(footer4);
			}


			function createFooter(footerParent,footerSrc,footerId){
					var footerData		=	footerSrc;
					var footerId			=	footerId;
					var footerBucketDiv		=	document.createElement('div');
						footerBucketDiv.id	=	footerId;		
					var footerBucketUl		=	document.createElement('ul');
					
						//for(each in footerData){
							footerData.each(function(value,index){
							var topLevel		=	footerData[index];
							var topfooterName		=	topLevel[0];
							var topfooterLink		=	topLevel[1];
							var footerItem		=	document.createElement('li');
							var footerItemA		=	document.createElement('a');
								footerItemA.href	=	topfooterLink;
								footerItemA.innerHTML	=	topfooterName;
							
							//footerp.appendChild(footerBucketDiv);
							footerBucketDiv.appendChild(footerBucketUl);
							footerBucketUl.appendChild(footerItem);
							footerItem.appendChild(footerItemA);
						});
						return footerBucketDiv;
				}

				function createLogo(logoSrc,logoClass){
						var logoData				=	logoSrc;
						var logoDisplay				=	document.createElement('a');
							logoDisplay.className	=	logoClass;
							logoDisplay.href		=	"/home/index.html";
						var logoSpan				=	document.createElement('span');
							logoSpan.className		=	"hide";		
							logoSpan.innerHTML		=	logoData;
								logoDisplay.appendChild(logoSpan);
								return logoDisplay;
					
				}




