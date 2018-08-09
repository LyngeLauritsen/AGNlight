function createForm()
{

var loginformDiv = $('login_form');

var loginformNoDrop = document.createElement('div');
loginformNoDrop.id="login_form_nodrop_old";

var loginlinks = document.createElement('span');
loginlinks.id = "login_links";

loginformDiv.appendChild(loginlinks);

/*var spanlogin = document.createElement('span');
spanlogin.innerHTML = "&rsaquo;&nbsp;";

var skipnavLogin = document.createElement('div');
skipnavLogin.className = "skiplinklogin";
skipnavLogin.innerHTML = '<a href="http://mynasa.nasa.gov/portal/site/mynasa/template.NASA_LOGIN_PROCESS">Follow this link to Login to MyNASA</a>';
var anchorlogin =  document.createElement('a');
anchorlogin.id = "loginnasa";
anchorlogin.className = "myOverlayLogin null bottom null observe_click";
anchorlogin.href = "#";
anchorlogin.innerHTML = "Log In To MyNASA";


spanlogin.appendChild(skipnavLogin);
spanlogin.appendChild(anchorlogin);


var textNode		=	document.createTextNode('|');


var spansingup = document.createElement('span');
spansingup.innerHTML = "&rsaquo;&nbsp;";

var anchorsignup = document.createElement('a');
anchorsignup.href="http://mynasa.nasa.gov/portal/site/mynasa/template.REGISTER";
anchorsignup.innerHTML = "Sign Up";

spansingup.appendChild(anchorsignup);


var ckUtil =  new CJL_CookieUtil("visitorinfo",0,"/",".nasa.gov");
var username=ckUtil.getSubValue("name");
var loginText		=	document.createElement('span');
loginText.innerHTML  = "Welcome "+username;



var logoutform = document.createElement('form'); 
logoutform.id = "gridLogout";
logoutform.name = "gridLogout";
logoutform.method = "post"
logoutform.action = "http://mynasa.nasa.gov/portal/site/mynasa/template.LOGOUT";


var spanlogout = document.createElement('span');
spanlogout.innerHTML = "&rsaquo;&nbsp;";


var logoutanchor = document.createElement('a');
logoutanchor.href = "javascript:gridLogoutSubmit();";
logoutanchor.innerHTML = "Log Out";

spanlogout.appendChild(logoutanchor);


var textNode1		=	document.createTextNode('|');
var textNode2		=	document.createTextNode('|');

var spanedit = document.createElement('span');
spanedit.innerHTML = "&rsaquo;&nbsp;";

var editanchor =document.createElement('a');
editanchor.href = "http://mynasa.nasa.gov/portal/site/mynasa/template.MY_ACCOUNT";
editanchor.innerHTML = "Edit Profile";

spanedit.appendChild(editanchor);;

var logouthidden = document.createElement('input');
logouthidden.type = "hidden";
logouthidden.id = "realm";
logouthidden.name = "realm";
logouthidden.value = "realml";



 if(ckUtil != null && username != null && username != ''){
		
		loginlinks.appendChild(loginText);
		loginlinks.appendChild(textNode1);
		loginlinks.appendChild(spanlogout);
		loginlinks.appendChild(textNode2);
		loginlinks.appendChild(spanedit);

		logoutform.appendChild(logouthidden);
		logoutform.appendChild(loginlinks);
		loginformDiv.innerHTML = "";
		loginformDiv.appendChild(logoutform);
	} else {
		loginlinks.appendChild(spanlogin);
		loginlinks.appendChild(textNode);
		loginlinks.appendChild(spansingup);

		loginformDiv.innerHTML = "";
		loginformDiv.appendChild(loginlinks) ;
	}*/
			

}

function gridLoginSubmit()
{
	var gridLoginform =$('gridLogin');
	gridLoginform.submit();
}

function gridLogoutSubmit()
{
	var gridLoginform =$('gridLogout');
	gridLogout.submit();
}

var text = false;
var textAllow = (window.location.search.indexOf('text=1')+1) ? false : true;

// User Preferences code ends

function switchText() {
  var val = (textAllow) ? '1' : '0';
	var s = window.location.href.split('#')[0];
  if(window.location.search) {
	  if(s.indexOf('text=')+1) {
		  s = s.replace('text='+s.split('text=')[1].split('&')[0],'text='+val);
		} else {
		  s += '&text='+val;
		}
	} else {
	  s += '?text='+val;
	}
	window.location.href = s;
}


function newAllowText() {
	//alert("Text");
  var s = '<span class="hide">&rsaquo;&nbsp;<a href="#" onclick="switchText(); return false;" >';
	s += (textAllow) ? 'Text Version' : 'Non-Text Version' ;
	s += '</a><br/><span>';
	return s;
}


function createFooterContent(editor,date,official,contact,link,sitemap) {
	if($('footer')!=null){
		if($('footercol1')){
			($('footercol1')).innerHTML = 'Page Last Updated:'+ date+' </br>'+
									'Page Editor:'+ editor+' <br />'+
									  'NASA Official:'+ official;
		}
		if($('footercol4')){
			  var liInnerHTML = ($('footercol4')).firstDescendant().firstDescendant();  /*This gives the contact li tag*/
			    var siteliTag    = ($('footercol4')).firstDescendant().firstDescendant().next(); /*This gives the sitemap li tag*/
			  var aInnerHTML = ($('footercol4')).firstDescendant().firstDescendant().firstDescendant();/*This gives the contact li anchor tag*/
			
			  var siteaTag    = ($('footercol4')).firstDescendant().firstDescendant().next().firstDescendant(); /*This gives the sitemap li anchor tag*/
			  
			  aInnerHTML.innerHTML=contact;
			  aInnerHTML.href=link;
             liInnerHTML.appendChild(aInnerHTML);

			 siteaTag.href=sitemap;
			 siteliTag.appendChild(siteaTag);
        }
	}
}



function CJL_CookieUtil(name, duration, path, domain, secure)
{
   this.affix = "";
   
   if( duration )
   {   	  
      var date = new Date();
	  var curTime = new Date().getTime();

	  date.setTime(curTime + (1000 * 60 * duration));
	  this.affix = "; expires=" + date.toGMTString();
   }
   
   if( path )
   {
      this.affix += "; path=" + path;
   }
   
   if( domain )
   {
      this.affix += "; domain=" + domain;
   }

   if( secure )
   {
      this.affix += "; secure=" + secure;
   }
   
      
   function getValue()
   {
      var m = document.cookie.match(new RegExp("(" + name + "=[^;]*)(;|$)"));

      return m ? m[1] : null;   
   }
   
   this.cookieExists = function()
   {
      return getValue() ? true : false;
   }
      
   this.expire = function()
   {
      var date = new Date();
	  date.setFullYear(date.getYear() - 1);
	  document.cookie=name + "=noop; expires=" + date.toGMTString(); 
   }
        
   this.setSubValue = function(key, value)
   {
      var ck = getValue();

      if( /[;, ]/.test(value) )
      {
         //Mac IE doesn't support encodeURI
		 value = window.encodeURI ? encodeURI(value) : escape(value);
      }

      
      if( value )
      {
         var attrPair = "@" + key + value;

         if( ck )
         {
             if( new RegExp("@" + key).test(ck) )
	         {
		        document.cookie =
				   ck.replace(new RegExp("@" + key + "[^@;]*"), attrPair) + this.affix;
	         }
	         else
	         {
		        document.cookie =
				   ck.replace(new RegExp("(" + name + "=[^;]*)(;|$)"), "$1" + attrPair) + this.affix;
	         }
         }
         else
         {
	        document.cookie = name + "=" + attrPair + this.affix;
         }
      }
      else
      {      
	     if( new RegExp("@" + key).test(ck) )
	     {
	        document.cookie = ck.replace(new RegExp("@" + key + "[^@;]*"), "") + this.affix;
	     }
      }
   }

      
   this.getSubValue = function(key)
   {
      var ck = getValue();

      if( ck )
      {
         var m = ck.match(new RegExp("@" + key + "([^@;]*)"));

	     if( m )
	     {
	        var value = m[1];

	        if( value )
	        { 
	           //Mac IE doesn't support decodeURI
			   return window.decodeURI ? decodeURI(value) : unescape(value);
	        }
	     }
      }
   }
}


function searchformsubmit() {
	var searchform = document.getElementById("search");
	
	if($("dropdown_search_label")!=null){
		var centername = $("dropdown_search_label").innerHTML.toLowerCase();
		if(centername=="nasa.gov"){
			document.getElementById("centername").value = "";
			searchform.action="http://search.nasa.gov/search/search.jsp";
		}else{
			document.getElementById("centername").value = centername;
			searchform.action="http://search.nasa.gov/search/centersearch.jsp?centername="+centername;
		}
	}else{
		searchform.action = "http://search.nasa.gov/search/search.jsp";
	}
	searchform.submit();
}

/*------------------New javascript for login and search -------------------- */
function createLoginForm() {

    var headerform = document.getElementById('header_form');

    var loginformDiv = document.createElement('div');
    loginformDiv.id = "login_form";

    var loginformDivNew = document.createElement('div');
    loginformDivNew.id = "login_form_new";

    var loginformNoDrop = document.createElement('div');
    loginformNoDrop.id = "login_form_nodrop";

	
	if (search_list.size() > 0) {
		headerform.appendChild(loginformDivNew);
	} else {
		headerform.appendChild(loginformNoDrop);
	}

    /*var loginlinks = document.createElement('span');
    loginlinks.id = "login_links";

    var spanlogin = document.createElement('span');
    spanlogin.innerHTML = "&rsaquo;&nbsp;";

    var skipnavLogin = document.createElement('div');
    skipnavLogin.className = "skiplinklogin";
    skipnavLogin.innerHTML = '<a href="http://mynasa.nasa.gov/portal/site/mynasa/template.NASA_LOGIN_PROCESS">Follow this link to Login to MyNASA</a>';
    var anchorlogin = document.createElement('a');
    anchorlogin.id = "loginnasa";
    anchorlogin.className = "myOverlayLogin null bottom null observe_click";
    anchorlogin.href = "#";
    anchorlogin.innerHTML = "Log In To MyNASA";

    spanlogin.appendChild(skipnavLogin);
    spanlogin.appendChild(anchorlogin);

    var textNode = document.createTextNode('|');

    var spansingup = document.createElement('span');
    spansingup.innerHTML = "&rsaquo;&nbsp;";

    var anchorsignup = document.createElement('a');
    anchorsignup.href = "http://mynasa.nasa.gov/portal/site/mynasa/template.REGISTER";
    anchorsignup.innerHTML = "Sign Up";

    spansingup.appendChild(anchorsignup);

    var ckUtil = new CJL_CookieUtil("visitorinfo", 0, "/", ".nasa.gov");
    var username = ckUtil.getSubValue("name");
    var loginText = document.createElement('span');
    loginText.innerHTML = "Welcome " + username;

    var logoutform = document.createElement('form');
    logoutform.id = "gridLogout";
    logoutform.name = "gridLogout";
    logoutform.method = "post"
    logoutform.action = "http://mynasa.nasa.gov/portal/site/mynasa/template.LOGOUT";

    var spanlogout = document.createElement('span');
    spanlogout.innerHTML = "&rsaquo;&nbsp;";

    var logoutanchor = document.createElement('a');
    logoutanchor.href = "javascript:gridLogoutSubmit();";
    logoutanchor.innerHTML = "Log Out";

    spanlogout.appendChild(logoutanchor);

    var textNode1 = document.createTextNode('|');
    var textNode2 = document.createTextNode('|');

    var spanedit = document.createElement('span');
    spanedit.innerHTML = "&rsaquo;&nbsp;";

    var editanchor = document.createElement('a');
    editanchor.href = "http://mynasa.nasa.gov/portal/site/mynasa/template.MY_ACCOUNT";
    editanchor.innerHTML = "Edit Profile";

    spanedit.appendChild(editanchor);;

    var logouthidden = document.createElement('input');
    logouthidden.type = "hidden";
    logouthidden.id = "realm";
    logouthidden.name = "realm";
    logouthidden.value = "realml";

    if (ckUtil != null && username != null && username != '') {
        loginlinks.appendChild(loginText);
        loginlinks.appendChild(textNode1);
        loginlinks.appendChild(spanlogout);
        loginlinks.appendChild(textNode2);
        loginlinks.appendChild(spanedit);

        logoutform.appendChild(logouthidden);
        logoutform.appendChild(loginlinks);
        if (search_list.size() > 0) {
            loginformDivNew.appendChild(logoutform)
            headerform.innerHTML = "";
            headerform.appendChild(loginformDivNew);
        } else {
            loginformDiv.appendChild(logoutform)
            headerform.innerHTML = "";
            headerform.appendChild(loginformDiv);
        }

    } else {
        loginlinks.appendChild(spanlogin);
        loginlinks.appendChild(textNode);
        loginlinks.appendChild(spansingup);

        if (search_list.size() > 0) {

            loginformDivNew.appendChild(loginlinks);
            headerform.appendChild(loginformDivNew);

        } else {
            loginformNoDrop.appendChild(loginlinks);
            headerform.appendChild(loginformNoDrop);

        }

    }*/

    if (search_list.size() > 0) {
        var searchSelect = new Element("select", {
            'disabled': "disabled"
        });
        search_list.each(function (searchList) {
            var opElem = new Element("option", {
                'id': searchList['id'],
                'name': searchList['name']
            });
            opElem.update(searchList['value']);
            searchSelect.appendChild(opElem);
        });
        if ($('searchselector') != null) {
            $('searchselector').appendChild(searchSelect);
            var skinnedDropper = new SkinnedSelectSearch($$('#searchselector')[0], $$('#searchselector' + ' select')[0], function () {}, '', 'gray');
        }
    }

}


function createSearchForm(){

	var headerform = document.getElementById('header_form');
	var searchformnasa = document.createElement('form');
	searchformnasa.id = "search";
	searchformnasa.method = "get";
	searchformnasa.action = "javascript:searchformsubmit();";

	var searchformcenter = document.createElement('form');
	searchformcenter.id = "search";
	searchformcenter.method = "get";
	searchformcenter.action = "javascript:searchformsubmit();";

	var searchdiv =  document.createElement('div');
	searchdiv.id = "search_form_new";

	var hiddenCenter = document.createElement('label');
	hiddenCenter.htmlFor = "searchfield";
	hiddenCenter.id = "searchfieldCenter";
	hiddenCenter.setAttribute('name','searchfieldCenter');
	hiddenCenter.innerHTML = '<input id="centername" name="centername" type="hidden" value=""/>';

	var searchdivNoDrop = document.createElement('div');
	searchdivNoDrop.id = "search_form_nodrop";
	
	var spansearchbtn = document.createElement('span');
	spansearchbtn.id = "searchbutton";

	var searchselector = document.createElement('div');
	searchselector.id = "searchselector";

	var spaninput = document.createElement('span');
	spaninput.id = "inputfield";

	var searchinput = document.createElement('input');
		searchinput.title = "searchfield";
		searchinput.type = "text";
		searchinput.id = "nasaInclude";
		searchinput.name = "nasaInclude";
		searchinput.className = "searchbox";
		searchinput.value = "";


	spaninput.appendChild(searchinput);

	var anchorbtn = document.createElement('a');
	anchorbtn.alt = "Search";
		anchorbtn.src = "/templateimages/redesign/modules/header/search-button.gif";
		anchorbtn.className = "searchbtn";
		anchorbtn.href="javascript:searchformsubmit();";


	spansearchbtn.appendChild(anchorbtn);

	var existingHeader = headerform.innerHTML;

	if(search_list.size()>0){
		searchdiv.appendChild(spaninput);
		searchdiv.appendChild(searchselector);
		searchdiv.appendChild(spansearchbtn);

		searchformcenter.appendChild(hiddenCenter);
		searchformcenter.appendChild(searchdiv);


		headerform.appendChild(searchformcenter);
	}
	else{
		searchdivNoDrop.appendChild(spaninput);
		searchdivNoDrop.appendChild(searchselector);
		searchdivNoDrop.appendChild(spansearchbtn);

		searchformnasa.appendChild(searchdivNoDrop);
		headerform.appendChild(searchformnasa);
	}



	


	if(search_list.size()>0){
	var searchSelect = new Element("select");
	search_list.each(function(searchList)
	{
		var opElem = new Element("option",{'id':searchList['id'],'name':searchList['name']});
		opElem.update(searchList['value']);
		searchSelect.appendChild(opElem);
	});
	if($('searchselector')!=null){
	$('searchselector').appendChild(searchSelect);
	var skinnedDropper = new SkinnedSelectSearch($$('#searchselector')[0],$$('#searchselector'+' select')[0],function(){},'','gray');}
	}
}