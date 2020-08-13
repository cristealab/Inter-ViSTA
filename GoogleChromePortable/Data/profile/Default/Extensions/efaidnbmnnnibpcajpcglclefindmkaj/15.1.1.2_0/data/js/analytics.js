var def;require=function(e){"use strict";return e},def=window.define?window.define:function(e,t){"use strict";return t.apply(null,[{ajax:$.ajax.bind($)}])};var exports=acom_analytics={};function dependOn(){"use strict";return[require("util"),require("proxy")]}function error(e){"use strict";SETTINGS.DEBUG_MODE&&console.error(e)}def(dependOn(),function(i,e){"use strict";var t,r=null,o=!0;for(t in r=r||new function(){e&&(this.proxy=e.proxy.bind(this)),this.e={EXTENSION_INSTALLED:"DCBrowserExt:Extension:Installed:Op",EXTENSION_UPDATE:"DCBrowserExt:Extension:Update:Op",EXTENSION_STARTUP:"DCBrowserExt:Extension:Startup:Op",EXTENSION_FORCE_UNINSTALL:"DCBrowserExt:Extension:ForceUninstall:Op",SIGN_IN_SHOWN:"DCBrowserExt:SignIn:OPERATION:Shown",SIGN_IN_COMPLETE:"DCBrowserExt:SignIn:OPERATION:Complete:Op",SIGN_IN_ABANDONED:"DCBrowserExt:SignIn:Abandoned:Op",SIGN_OUT_CLICKED:"DCBrowserExt:SignOut:Clicked",FLICKR_OFFER_SHOWN:"DCBrowserExt:Flickr:Offer:Shown",FLICKR_OFFER_CLICKED:"DCBrowserExt:Flickr:Offer:Clicked",FLICKR_CONTEXT_CLICK:"DCBrowserExt:Flickr:Context:Clicked",CONTEXT_UPLOAD_PDF_PAGE:"DCBrowserExt:Context:Upload:PdfPage:Clicked",CONTEXT_UPLOAD_LINK:"DCBrowserExt:Context:Upload:Link:Clicked",CONTEXT_UPLOAD_IMAGE:"DCBrowserExt:Context:Upload:Image:Clicked",CONTEXT_MENU_CONVERT_PAGE:"DCBrowserExt:ContextMenu:ConvertPage:Clicked",CONTEXT_MENU_CONVERT_LINK:"DCBrowserExt:ContextMenu:ConvertLink:Clicked",CONTEXT_MENU_APPEND_PAGE:"DCBrowserExt:ContextMenu:AppendPage:Clicked",CONTEXT_MENU_APPEND_LINK:"DCBrowserExt:ContextMenu:AppendLink:Clicked",REDIRECT:"DCBrowserExt:Redirect:OPERATION:Op",PDF_MENU_UPLOAD_COMPLETE_SHOWN:"DCBrowserExt:PDF:Menu:OPERATION:Upload:Complete:Shown",PDF_MENU_UPLOAD_CLICKED:"DCBrowserExt:PDF:Menu:Copy:Clicked",PDF_MENU_EXPORT_CLICKED:"DCBrowserExt:PDF:Menu:Export:Clicked",PDF_MENU_SEND_CLICKED:"DCBrowserExt:PDF:Menu:Send:Clicked",PDF_MENU_FILLSIGN_CLICKED:"DCBrowserExt:PDF:Menu:FillSign:Clicked",PDF_SOURCE_SIZE:"DCBrowserExt:PDF:Source:Size:RANGE:Op",TREFOIL_CLICKED:"DCBrowserExt:Trefoil:Clicked",TREFOIL_PDF_MENU_SHOWN:"DCBrowserExt:TrefoilMenu:PDF:Popup:Shown",TREFOIL_PDF_FROM_CLICK:"DCBrowserExt:TrefoilMenu:PDF:FromClick:TIREKICK:Shown",TREFOIL_PDF_VISIT_AIC:"DCBrowserExt:TrefoilMenu:PDF:VisitAIC:TIREKICK:Clicked",TREFOIL_PDF_ACROBAT:"DCBrowserExt:TrefoilMenu:PDF:OpenInAcrobat:TIREKICK:Clicked",TREFOIL_PDF_READER:"DCBrowserExt:TrefoilMenu:PDF:OpenInReader:TIREKICK:Clicked",TREFOIL_PDF_DOWNLOAD_OPENED:"DCBrowserExt:TrefoilMenu:PDF:OpenInAcrobat:Complete:Op",TREFOIL_PDF_DOWNLOAD_OPENED_READER:"DCBrowserExt:TrefoilMenu:PDF:OpenInReader:Complete:Op",TREFOIL_PDF_DOWNLOAD_OPEN_FAILED:"DCBrowserExt:TrefoilMenu:PDF:OpenInAcrobat:Failed:Op",TREFOIL_PDF_DOWNLOAD_OPEN_FAILED_READER:"DCBrowserExt:TrefoilMenu:PDF:OpenInReader:Failed:Op",TREFOIL_HTML_FROM_CLICK:"DCBrowserExt:TrefoilMenu:HTML:FromClick:Shown",TREFOIL_HTML_OPTIONS_FROM_CLICK:"DCBrowserExt:TrefoilMenu:HTML:WithOptions:FromClick:Shown",TREFOIL_HTML_VISIT_AIC:"DCBrowserExt:TrefoilMenu:HTML:VisitAIC:TIREKICK:Clicked",TREFOIL_HTML_PREFERENCES_CLICK:"DCBrowserExt:TrefoilMenu:HTML:Preferences:Clicked",TREFOIL_HTML_CONVERT_NEW:"DCBrowserExt:TrefoilMenu:HTML:ConvertHTML:New:TIREKICK:Clicked",TREFOIL_HTML_CONVERT_APPEND:"DCBrowserExt:TrefoilMenu:HTML:ConvertHTML:Append:TIREKICK:Clicked",TREFOIL_HTML_CONVERT_NO_OPEN:"DCBrowserExt:TrefoilMenu:HTML:ConvertHTML:NoOpen:Clicked",TREFOIL_HTML_CONVERT_OPEN_DEFAULT:"DCBrowserExt:TrefoilMenu:HTML:ConvertHTML:WithOpen:Default:Clicked",TREFOIL_HTML_CONVERT_OPEN_CHANGED:"DCBrowserExt:TrefoilMenu:HTML:ConvertHTML:WithOpen:Changed:Clicked",TREFOIL_HTML_CONVERT_WAITING:"DCBrowserExt:TrefoilMenu:HTML:ConvertHTML:Waiting:Shown",TREFOIL_HTML_CONVERT_DOWNLOADING:"DCBrowserExt:TrefoilMenu:HTML:ConvertHTML:Downloading:Shown",TREFOIL_HTML_CONVERT_IN_PROGRESS:"DCBrowserExt:TrefoilMenu:HTML:ConvertHTML:InProgress:Shown",TREFOIL_HTML_CONVERT_CANCELLED:"DCBrowserExt:TrefoilMenu:HTML:ConvertHTML:Cancelled:Shown",TREFOIL_HTML_CONVERT_COMPLETE:"DCBrowserExt:TrefoilMenu:HTML:ConvertHTML:Complete:Shown",TREFOIL_HTML_CONVERT_FAILED:"DCBrowserExt:TrefoilMenu:HTML:ConvertHTML:Failed:Shown",TREFOIL_HTML_CONVERT_NO_ACROBAT:"DCBrowserExt:TrefoilMenu:HTML:ConvertHTML:NoAcrobat:Shown",TREFOIL_HTML_OPENPDF_PREF_ON:"DCBrowserExt:TrefoilMenu:HTML:OpenPDFPref:On:Clicked",TREFOIL_HTML_OPENPDF_PREF_OFF:"DCBrowserExt:TrefoilMenu:HTML:OpenPDFPref:Off:Clicked",PERSIST_PDF_MENU_SHOWN:"DCBrowserExt:PersistMenu:PDF:Popup:Shown",PERSIST_PDF_ACROBAT:"DCBrowserExt:PersistMenu:PDF:OpenInAcrobat:TIREKICK:Clicked",PERSIST_PDF_DOWNLOAD_OPENED:"DCBrowserExt:PersistMenu:PDF:OpenInAcrobat:Complete:Op",PERSIST_PDF_DOWNLOAD_OPEN_FAILED:"DCBrowserExt:PersistMenu:PDF:OpenInAcrobat:Failed:Op",PERSIST_PDF_MENU_CLOSED:"DCBrowserExt:PersistMenu:PDF:Close:TIREKICK:Clicked",PERSIST_PDF_OPENPDF_PREF_OFF:"DCBrowserExt:PersistMenu:PDF:OpenPDFPref:Off:Clicked",UPLOAD_PROGRESS_SHOWN:"DCBrowserExt:Upload:OPERATION:Progress:Shown",CREATE_FORM_PROGRESS_SHOWN:"DCBrowserExt:Upload:FillSign:CreateForm:Progress:Shown",UPLOAD_COMPLETE:"DCBrowserExt:Upload:OPERATION:Complete:Op",CREATE_FORM_COMPLETE:"DCBrowserExt:Upload:FillSign:CreateForm:Complete:Op",UPLOAD_RENAME_CLICKED:"DCBrowserExt:Upload:RenameOrMove:Clicked",ERROR_SHOWN:"DCBrowserExt:Error:Shown",ERROR_WRONG_MIME_TYPE:"DCBrowserExt:Error:WrongMimeType:Shown",OPTIONS_SET_ENV:"DCBrowserExt:Options:SetEnv:ENVIRONMENT:Op",OPTIONS_ENABLE_HTML2PDF:"DCBrowserExt:Options:EnableHTML2PDF:ENVIRONMENT:Op",HTML_SOURCE_SIZE:"DCBrowserExt:HTML:Source:Size:RANGE:Op",HTML_SOURCE_SIZE_TOO_LARGE_ERROR:"DCBrowserExt:HTML:Source:Size:TooLarge:Error:Shown",HTML_SOURCE_SITE:"DCBrowserExt:HTML:Source:Site:SITE",HTML_SOURCE_PROTOCOL:"DCBrowserExt:HTML:Source:Protocol:PROTOCOL",HTML_SOURCE_CONTENT:"DCBrowserExt:HTML:Source:CONTENT:Op",HTML_CONVERSION_STAGE_TIMING:"DCBrowserExt:HTML:Conversion:STAGE:TIMING",OS_MAC_OP:"DCBrowserExt:OS:mac:Op",OS_WIN_OP:"DCBrowserExt:OS:win:Op",SHIM_VERSION:"DCBrowserExt:Shim:Version:VERSION:Op",OPTIONS_PAGE:"DCBrowserExt:OptionsPage:Shown",FTE_LAUNCH:"DCBrowserExt:FTE:Launch:Shown"},this.event=function(e,t){if(SETTINGS.ANALYTICS_OPT_IN&&o)if(e){t=t||{};try{this.s.t(function(e,r,t){return r=r.replace(/OPERATION/,e.operation||"Unknown").replace(/ENVIRONMENT/,e.environment).replace(/STAGE/,e.stage).replace(/TIMING/,e.timing).replace(/RANGE/,e.size).replace(/SITE/,e.site).replace(/PROTOCOL/,e.protocol),i.each(t,function(e,t){r=r.replace(e,t)}),SETTINGS.ANALYTICS&&console.log("%c"+r,"color: #800080"),{eVar1:e.version,eVar2:e.installType,eVar3:e.environment,eVar4:e.shim,pageName:r}}(this,e,t))}catch(e){error(e.toString())}}else error("Missing analytics string")},this.init=function(e,t){var r=i.getCookie("logAnalytics");o=!r||"true"===i.getCookie("logAnalytics"),this.shim="not_set",this.version=e,this.installType=t,this.environment="prod",this.s=new AppMeasurement,this.s.ssl=!0,"development"===this.installType?this.s.account="adbcreatepdfplugin.dev":this.s.account="adbcreatepdfplugin.prod",this.s.trackingServer="stats.adobe.com",this.s.trackingServerSecure="sstats.adobe.com"},this.setArg=function(e,t){e&&(this[e]=t)},this.setParamsAndLogAnalytics=function(e,t,r){if(e){var o=this;e.forEach(function(e){o[r]=e,o.event(t)})}},this.setAnalyticsUsage=function(e,t){o=e,i.setCookie("logAnalytics",o.toString(),3650),i.sendMessage({options_op:"saved_analytics",tabId:t})},this.getAnalyticsUsage=function(){return o},this.setOp=function(e){e&&(this.operation=e)},this.error=function(e,t,r){try{this.event(this.e.ERROR_SHOWN);var o="DCBrowser:Error:JS:"+r.stack.match(/([A-Za-z0-9\-]+)\.js:(\d*):(\d*)/)[0]+":"+r.message.replace(/ /g,"_");this.event(o)}catch(e){}},e&&e.handlers(this.error.bind(this)),this.checkSizes=function(e){0<e&&e<=1?this.setArg("size","0_1"):1<e&&e<=2?this.setArg("size","1_2"):2<e&&e<=5?this.setArg("size","2_5"):5<e&&e<=10?this.setArg("size","5_10"):10<e&&e<=50?this.setArg("size","10_50"):50<e&&e<=500?this.setArg("size","50_500"):500<e&&e<=1e3?this.setArg("size","500_1000"):1e3<e&&e<=2e3?this.setArg("size","1000_2000"):2e3<e&&e<=3e3?this.setArg("size","2000_3000"):3e3<e&&e<=4e3?this.setArg("size","3000_4000"):4e3<e&&this.setArg("size","4000_")},this.checkAndLogHTMLBlobSize=function(e){this.checkSizes(e),this.event(this.e.HTML_SOURCE_SIZE)},this.checkAndLogPDFSize=function(e){this.checkSizes(e),this.event(this.e.PDF_SOURCE_SIZE)},this.checkAndLogTimingRange=function(e){0<e&&e<=5?this.setArg("timing","0_5"):5<e&&e<=10?this.setArg("timing","5_10"):10<e&&e<=20?this.setArg("timing","10_20"):20<e&&e<=30?this.setArg("timing","20_30"):30<e&&e<=50?this.setArg("timing","30_50"):50<e&&e<=100?this.setArg("timing","50_100"):100<e&&e<=200?this.setArg("timing","100_200"):200<e&&e<=600?this.setArg("timing","200_600"):600<e&&e<=1200?this.setArg("timing","600_1200"):1200<e&&e<=3e3?this.setArg("timing","1200_3000"):3e3<e&&this.setArg("timing","3000_"),this.event(this.e.HTML_CONVERSION_STAGE_TIMING)},this.logSiteAndProtocolAnalytics=function(e){var t,r,o,i,n=[],s=[];0!==e.indexOf("chrome:")&&(o=this.parseURL(e),r=/^(http|https):\/\/www/.test(e),t=/^(([0-9]+\.){3}([0-9]+))$/.test(e),i=o.hostname.split(".").reverse()[0].toLowerCase(),r&&n.push("WWWW"),t&&n.push("IP"),"com"!==i&&"org"!==i&&"net"!==i&&"int"!==i&&"edu"!==i&&"gov"!==i&&"mil"!==i||n.push("TLD"),s.push(o.protocol),this.setParamsAndLogAnalytics(s,this.e.HTML_SOURCE_PROTOCOL,"protocol"),this.setParamsAndLogAnalytics(n,this.e.HTML_SOURCE_SITE,"site"))},this.parseURL=function(e){var t=document.createElement("a");return t.href=e,{protocol:t.protocol.slice(0,t.protocol.length-1),host:t.host,hostname:t.hostname,port:t.port,pathname:t.pathname}},this.logBrowserAnalytics=function(e){e.analytics&&(e.analytics.forEach(this.proxy(function(e){this.checkAndLogAnalytics(e)})),delete e.analytics)},this.logContents=function(e){e.content_analytics&&(e.content_analytics.forEach(this.proxy(function(e){this.event(this.e.HTML_SOURCE_CONTENT,{CONTENT:e})})),delete e.analytics)},this.checkAndLogAnalytics=function(e){var t,r="TIREKICK",o={};e?(-1!==e.indexOf(r)&&((t=i.getCookie(r))?-1!==t.indexOf(e)?o[r]="Subsequent":(i.setCookie(r,t+"|"+e,3650),o[r]="FirstTime"):(i.setCookie(r,e,3650),o[r]="FirstTime")),this.event(e,o)):error("Missing analytics string")},this.logError=function(e){var t=e;"web2pdfHTMLTooLarge"===e&&(t=this.e.HTML_SOURCE_SIZE_TOO_LARGE_ERROR),this.event(t)}})r.hasOwnProperty(t)&&("function"==typeof r[t]?exports[t]=r[t].bind(r):exports[t]=r[t]);return r});