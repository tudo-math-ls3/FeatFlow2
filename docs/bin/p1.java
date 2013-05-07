import java.io.*;
import java.util.*;
import java.util.regex.Pattern;

import javax.xml.parsers.*;
import javax.xml.transform.*;
import javax.xml.transform.dom.DOMSource;
import javax.xml.transform.stream.StreamResult;

import org.w3c.dom.*;
import org.xml.sax.*;

import java.net.URLEncoder;

class p1
{

  static String progname = "p1.class";

  static FileOutputStream writeFile;
  static DataOutputStream out;

  static int modus;
  static String filename;

  static int objidx = 0;
  static int objentidx = 1;
  static int attridx = 1;

  static final int FORMAT_HTML = 0;
  static final int FORMAT_TEX = 1;
  static final int FORMAT_TEXMOD = 2;

  public static void parse_vars(String value2, int mode)
  {
    try
    {

      String collect = "";

      String[] lines = value2.split("\n");

      int tmodus;
      String filler;

      if (mode==0)
      {
        tmodus=1;
        filler="";
      }
      else
      {
        filler="";
        switch(modus)
        {
          case FORMAT_HTML:
            filler="<td>&nbsp;</td>";
            break;
          case FORMAT_TEX:
            filler="&";
            break;
        }

        tmodus=0;
      }


      for(int l = 0; l<lines.length; l++)
      {

        if (lines[l].trim().indexOf("!")==0)
        {
          collect=collect.concat(" "+lines[l].trim().substring(1));
        }
        else
        {
          if ( (tmodus==0) && (lines[l].toLowerCase().indexOf("type")>0) && 
             ( (lines[l].toLowerCase().indexOf("type (")==-1) && 
             (lines[l].toLowerCase().indexOf("type(")==-1)))
          {
            int cpol1=lines[l].toLowerCase().indexOf("type");

            String ident = lines[l].substring(cpol1+4);

            tmodus=1;

            switch(modus)
            {
              case FORMAT_HTML:
                out.writeBytes("<tr><td><a name="+ident+"><em>"+ident+"</em></a></td><td colspan=4>"+collect+"</td></tr>\n");
                break;
              case FORMAT_TEX:
                out.writeBytes("{\\tt "+ident.replaceAll("character","char").replaceAll("_","\\\\_").replaceAll(" ->"," \\\\ensuremath{\\\\rightarrow}")
                              +"} & {"+collect.replaceAll("_","\\\\_").replaceAll(" ->"," \\\\ensuremath{\\\\rightarrow}")+"}\\\\\n");
                break;
            }

            collect="";
          }

          if ( (tmodus==1) && (lines[l].toLowerCase().indexOf("end type")>0))
          {
            tmodus=0;

            switch(modus)
            {
              case FORMAT_HTML:
                out.writeBytes("<tr><td colspan=5>&nbsp;</td><td>\n");
                break;
            }
          }

          if (tmodus==1)
          {

            int cpol1=lines[l].indexOf("::");

            if (cpol1>0)
            {
              String typefull = lines[l].substring(1,cpol1).trim();

              String[] attr = typefull.split(",");

              String vtype=attr[0];
              String rank = "";
              String name;

              if (attr.length>1)
              {
                if (attr[1].toLowerCase().indexOf("dimension")>0)
                {
                  rank = attr[1];
                  for (int ia=2;ia<attr.length;ia++)
                  {
                    if ((attr[ia].indexOf(")")>0) && (attr[ia].toLowerCase().indexOf("intent")==-1))
                      rank = rank.concat(", "+attr[ia]);
                    else if (attr[ia].matches("0-9"))
                      rank = rank.concat(", "+attr[ia]);
                  }
                }
                else
                {
                  rank="";
                  switch(modus)
                  {
                    case FORMAT_HTML:
                      rank = "&nbsp;";
                      break;
                  }
                }
              }
              else
              {
                rank="";

                switch(modus)
                {
                  case FORMAT_HTML:
                    rank = "&nbsp;";
                    break;
                }
              }

              name=lines[l].substring(cpol1+2).trim();

              if (name.trim().charAt(0)=='c')
              {
                switch(modus)
                {
                  case FORMAT_HTML:
                    out.writeBytes("<tr>"+filler+"<td><a name="+name+"></a><a href=\"#"+name+"linkdef\">"+name+"</a></td><td><em>"+vtype+"</em></td><td>"+rank+"</td><td>"+collect+"</td></tr>\n");
                    collect="";
                    break;
                  case FORMAT_TEX:
                    out.writeBytes(filler+" {\\tt "+name.replaceAll("_","\\\\_").replaceAll(" ->"," \\\\ensuremath{\\\\rightarrow}")
                    +"} & "+vtype.replaceAll("_","\\\\_").replaceAll(" ->"," \\\\ensuremath{\\\\rightarrow}").replaceAll("character","char")+
                      " & "+rank.replaceAll("_","\\\\_").replaceAll(" ->"," \\\\ensuremath{\\\\rightarrow}").replaceAll("dimension","dim")
                      +" & "+collect.replaceAll("_","\\\\_").replaceAll(" ->"," \\\\ensuremath{\\\\rightarrow}")+"\\\\ \\hline \n");
                    collect="";
                    break;
                }
              }
              else if (typefull.toLowerCase().indexOf("type (")>-1)
              {
                String kname=typefull.substring(typefull.indexOf("(")+1,typefull.indexOf(")")).trim();
                switch(modus)
                {
                  case FORMAT_HTML:
                    out.writeBytes("<tr>"+filler+"<td>"+name+"</td><td><a href=\"#"+kname+"\">"+kname+"</a></td><td>"+rank+"</td><td>"+collect+"</td></tr>\n");
                    collect="";
                    break;
                  case FORMAT_TEX:
                    out.writeBytes(filler+" {\\tt "+name.replaceAll("_","\\\\_").replaceAll(" ->"," \\\\ensuremath{\\\\rightarrow}")+"} & {\\tt "+
                    kname.replaceAll("_","\\\\_").replaceAll(" ->"," \\\\ensuremath{\\\\rightarrow}").replaceAll("character","char")
                      +"} & "+rank.replaceAll("_","\\\\_").replaceAll(" ->"," \\\\ensuremath{\\\\rightarrow}").replaceAll("dimension","dim")
                      +" &  "+collect.replaceAll("_","\\\\_").replaceAll(" ->"," \\\\ensuremath{\\\\rightarrow}")+"\\\\ \\hline \n");
                    collect="";
                    break;
                }
              }
              else
              {
                if (name.indexOf("=")>-1) name=name.substring(0,name.indexOf("="));

                switch(modus)
                {

                  case FORMAT_HTML:
                    out.writeBytes("<tr>"+filler+"<td>"+name+"</td><td><em>"+vtype+"</em></td><td>"+rank+"</td><td>"+collect+"</td></tr>\n");
                    collect="";
                    break;
                  case FORMAT_TEX:
                    out.writeBytes(filler+" {\\tt "+name.replaceAll("_","\\\\_").replaceAll(" ->"," \\\\ensuremath{\\\\rightarrow}")+"} & {\\tt "
                      +vtype.replaceAll("_","\\\\_").replaceAll(" ->"," \\\\ensuremath{\\\\rightarrow}").replaceAll("character","c").replaceAll("len=","")
                      +"} & "+rank.replaceAll("_","\\\\_").replaceAll(" ->"," \\\\ensuremath{\\\\rightarrow}").replaceAll("dimension","").replaceAll("(\\(|\\))","")
                      +" & "+collect.replaceAll("_","\\\\_").replaceAll(" ->"," \\\\ensuremath{\\\\rightarrow}")+"\\\\ \\hline\n");
                    collect="";
                    break;
                }
                collect="";
              }

            }
          }
        }
      }
    }
    catch (Exception ex)
    {
      System.out.println("Exp "+ex);
      ex.printStackTrace();
    }
  }

  public static void main(String args[])
  {
    Pattern patternFunction   = Pattern.compile(".*\\bfunction\\b.*", Pattern.CASE_INSENSITIVE | Pattern.DOTALL);
    Pattern patternSubroutine = Pattern.compile(".*\\bsubroutine\\b.*", Pattern.CASE_INSENSITIVE | Pattern.DOTALL);

    // Determine parsing mode and name of file to parse
    switch (args.length)
    {
      // nothing given, default values for everything
      case 0:
	  System.out.println(progname + ": No output mode given, HTML output selected.");
	  modus = FORMAT_HTML;

	  System.out.println(progname + ": No input file given, <module.f90> assumed.");
	  filename = "module.f90";
	  break;

      // only parsing mode given, use default file name
      case 1:
	  if (args[0].equals("html"))
	  {
	    System.out.println(progname + ": HTML output selected.");
	    modus = FORMAT_HTML;
	  }
	  else if (args[0].equals("texmod"))
	  {
	    System.out.println(progname + ": TEX output (module description only) selected.");
	    modus = FORMAT_TEXMOD;
	  }
	  else if (args[0].equals("tex"))
          {
	    System.out.println(progname + ": TEX output selected.");
	    modus = FORMAT_TEX;
          }
	  else
	  {
	    System.out.println(progname + ": Unrecognised output format, HTML output selected.");
	    modus = FORMAT_HTML;
          }

	  filename = "module.f90";
	  System.out.println(progname + ": No input file given, <" + filename + "> assumed.");
	  break;

      // both given: parsing mode given and file name
      default:
	  if (args[0].equals("html"))
	  {
	    System.out.println(progname + ": HTML output selected.");
	    modus = FORMAT_HTML;
	  }
	  else if (args[0].equals("texmod"))
	  {
	    System.out.println(progname + ": TEX output (module description only) selected.");
	    modus = FORMAT_TEXMOD;
	  }
	  else if (args[0].equals("tex"))
          {
	    System.out.println(progname + ": TEX output selected.");
	    modus = FORMAT_TEX;
          }
	  else
	  {
	    System.out.println(progname + ": Unrecognised output format, HTML output selected.");
	    modus = FORMAT_HTML;
          }

	  filename = args[1];
	  System.out.println(progname + ": Using input file <" + filename + ">.");
	  break;
    }

    try
    {
      // Strip off extension from filename
      // (to be able to create an accordingly named output file)
      String basename = filename.substring(0, filename.lastIndexOf("."));

      File f = new File(filename);

      switch(modus)
      {
        case FORMAT_HTML:
          writeFile = new FileOutputStream(basename + ".html");
          out = new DataOutputStream(writeFile);
          break;
        case FORMAT_TEX:
          writeFile = new FileOutputStream(basename + ".docu.tex");
          out = new DataOutputStream(writeFile);
          break;
        case FORMAT_TEXMOD:
          writeFile = new FileOutputStream(basename + ".descr.tex");
          out = new DataOutputStream(writeFile);
          break;
      }

      // prepare XML parsing
      DocumentBuilderFactory factory =
        DocumentBuilderFactory.newInstance();
     // factory.setValidating(false);
     // factory.setNamespaceAware(true);
      factory.setIgnoringElementContentWhitespace(true);
      factory.setIgnoringComments(true);

      DocumentBuilder builder = factory.newDocumentBuilder();

      builder.setErrorHandler(new org.xml.sax.ErrorHandler()
      {
        // ignore fatal errors (an exception is guaranteed)
        public void fatalError(SAXParseException err)
          throws SAXException
        {
         System.out.println("** FATALERROR"
            + ", line "
            + err.getLineNumber()
            + ", uri "
            + err.getSystemId());


        }

        // treat validation errors as fatal
        public void error(SAXParseException err)
          throws SAXParseException
        {
         System.out.println("** ERROR"
            + ", line "
            + err.getLineNumber()
            + ", uri "
            + err.getSystemId());
        }

        // dump warnings too
        public void warning(SAXParseException err)
          throws SAXParseException
        {
          System.out.println("** Warning"
            + ", line "
            + err.getLineNumber()
            + ", uri "
            + err.getSystemId());
        }
      });

      // parse actual file
      Document doc = builder.parse("file:" + f.getAbsolutePath());

      // get root

      Element root = doc.getDocumentElement();
      // get all optionLists

      switch(modus)
      {
        case FORMAT_HTML:

//        out.writeBytes("<!doctype HTML PUBLIC \"-//W3C//DTD HTML 4.0 Transitional//EN\">\n");
//        out.writeBytes("<HTML>\n");
//        out.writeBytes("<HEAD>\n");
//        out.writeBytes("<META http-equiv=\"Content-Type\" content=\"text/html; charset=iso-8859-1\">\n");
//        out.writeBytes("<TITLE>FBDB Dokumentation</TITLE>\n");
//        out.writeBytes("<TITLE>FBDB Dokumentation</TITLE>\n");
//        out.writeBytes("<STYLE TYPE=\"text/css\" media=\"screen\">\n");
//        out.writeBytes("<!-- @import \"feast.css\"; -->\n");
//        out.writeBytes("</STYLE>\n");
//        out.writeBytes("</HEAD>\n");

//           out.writeBytes("<BODY>\n");
//           out.writeBytes("<H1>FEAT2 Dokumentation</H1>\n<br><br>");
          break;
      }

      NodeList lists = root.getChildNodes();
      NamedNodeMap attrg = root.getAttributes();

      int firstsr=0;
      int firstf=0;

      // iterate over lists and parse each
      for (int i = 0; i < lists.getLength(); i++)
      {

        Node element =  lists.item(i);

        String type = element.getNodeName();

        if (type.equals("globals"))
        {

          switch(modus)
          {
            case FORMAT_HTML:
              out.writeBytes("<br><hr><br><a name=\"GLOBAL\"> Global variable definitions </a> <br><br>\n");
              out.writeBytes("<table border=0>\n");
              out.writeBytes("<th>Name</th><th>Type</th><th>Rank</th><th>Purpose</th>\n");

              break;
            case FORMAT_TEX:
              out.writeBytes("{\\bf Global variable definitions: }\n");
              out.writeBytes("\\begin{center}\\small\n");
              out.writeBytes("\\phantom{fucksupertabular}\\begin{supertabular}{l l l p{5cm}}\n");
              out.writeBytes("Name&Type&Rank&Purpose\\\\ \\hline \n");
              break;

          }

          NodeList lists2 = element.getChildNodes();

          Node element2 =  lists2.item(0);

          String type2 = element2.getNodeName();

          String value2=element2.getNodeValue();

          parse_vars(value2,0);

          switch(modus)
          {
            case FORMAT_HTML:
              out.writeBytes("</table>\n");
              out.writeBytes("<hr>\n");
              break;
            case FORMAT_TEX:
              out.writeBytes("\\end{supertabular}\n");
              out.writeBytes("\\end{center}\n\n");
              break;

          }
        }

        //===========================================================================

        if (type.equals("name"))
        {
          NodeList lists1 = element.getChildNodes();

          String value1=lists1.item(0).getNodeValue().trim();
          switch(modus)
          {
            case FORMAT_HTML:
              out.writeBytes("<h1>Module "+value1+"</h1>\n");
              break;
            case FORMAT_TEX:
              out.writeBytes("\\section{Module "+value1+"}\n");
              break;
          }
        }

        if (type.equals("version"))
        {
          NodeList lists1 = element.getChildNodes();

          String value1=lists1.item(0).getNodeValue().trim();

          switch(modus)
          {
            case FORMAT_HTML:
              out.writeBytes("<h2>Version "+value1+"</h2>\n");
              break;
            case FORMAT_TEX:
              out.writeBytes("{\\bf Version: }"+value1+"\\\\[0.5cm]\n");
              break;
          }
        }

        if (type.equals("purpose"))
        {
          NodeList lists1 = element.getChildNodes();

          String value1 = lists1.item(0).getNodeValue();
	  value1 = value1.replaceAll("\n","");
	  value1 = value1.replaceAll("!#","");
	  value1 = value1.replaceAll("#","");
	  value1 = value1.trim();

          switch(modus)
          {
            case FORMAT_HTML:
              out.writeBytes("\n<div class=\"purpose\">\n    <h4>Purpose:</h4>\n   " + value1 + "</div>\n\n");
              break;
            case FORMAT_TEXMOD:
              out.writeBytes(value1.replaceAll("_","\\\\_").replaceAll(" ->"," \\\\ensuremath{\\\\rightarrow}") + "\n");
              break;
            case FORMAT_TEX:
              out.writeBytes("{\\bf Purpose:} " + value1.replaceAll("_","\\\\_").replaceAll(" ->"," \\\\ensuremath{\\\\rightarrow}") + "\\\\[0.5cm]\n");
              break;
          }
        }

        //===========================================================================


        if (type.equals("subroutine"))
        {
          if (firstsr==0)
          {
            firstsr=1;
            switch(modus)
            {
              case FORMAT_HTML:
                out.writeBytes("<hr>\n\n<div class=\"los\"><a name=\"SUBROUTINE\">List of subroutines</a></div>\n\n");
                break;
            }
          }

          NodeList lists1 = element.getChildNodes();
          for (int i1 = 0; i1 < lists1.getLength(); i1++)
          {
            Node element1 =  lists1.item(i1);
            String type1 = element1.getNodeName();

            if (type1.equals("#text"))
            {
              String value1=element1.getNodeValue();

	      if (patternSubroutine.matcher(value1).matches())
              {
                String mame = value1.substring(value1.toLowerCase().indexOf("subroutine")+10,value1.indexOf("(")).trim();

                switch(modus)
                {
                  case FORMAT_HTML:
                    out.writeBytes("<a href=\"#"+mame+"\">"+mame+"</a><br>\n");
                    break;
                }
              }
            }
          }
        }

        //===========================================================================

        if (type.equals("function"))
        {
          if (firstf==0)
          {
            firstf=1;
            switch(modus)
            {
              case FORMAT_HTML:
                out.writeBytes("<hr>\n\n<div class=\"lof\"><a name=\"FUNCTION\">List of functions</a></div>\n\n");
                break;
            }
          }

          NodeList lists1 = element.getChildNodes();
          for (int i1 = 0; i1 < lists1.getLength(); i1++)
          {
            Node element1 =  lists1.item(i1);
            String type1 = element1.getNodeName();

            if (type1.equals("#text"))
            {
              String value1=element1.getNodeValue();
	      if (patternFunction.matcher(value1).matches())
              {
                String mame = value1.substring(value1.toLowerCase().indexOf("function")+9,value1.lastIndexOf("(")).trim();

                switch(modus)
                {
                  case FORMAT_HTML:
                    out.writeBytes("<a href=\"#"+mame+"\">"+mame+"</a><br>\n");
                    break;
                }
              }
            }
          }
        }

        //===========================================================================


        if (type.equals("types"))
        {
          switch(modus)
          {
            case FORMAT_HTML:
              out.writeBytes("<br><hr><br><a name=\"TYPE\"> Type definitions </a> <br> <br>\n");
              out.writeBytes("<table border=0>\n");
              out.writeBytes("<th>Type name</th><th>Component name</th><th>Type</th><th width=128>Rank</th><th>Purpose</th>\n");
              break;
            case FORMAT_TEX:
              out.writeBytes("{\\bf Type definitions: }\n");
              out.writeBytes("\\begin{center}\\small\n");
              out.writeBytes("\\phantom{fucksupertabular}\\begin{supertabular}{l l l l p{3cm}}\n");
              out.writeBytes("Type name&component name&Type&Rank&Purpose\\\\ \\hline \n");
              break;
          }

          NodeList lists1 = element.getChildNodes();

          for (int i1 = 0; i1 < lists1.getLength(); i1++)
          {
            Node element1 =  lists1.item(i1);
            String type1 = element1.getNodeName();


            if (type1.equals("typeblock"))
            {

              switch(modus)
              {
                case FORMAT_HTML:
                  out.writeBytes("<tr><td colspan=5></td></tr>\n");
                  break;
                case FORMAT_TEX:
                 out.writeBytes("\\\\\n");
                 break;
              }

              NodeList lists2 = element1.getChildNodes();

              int modus=0;

              for (int i2 = 0; i2 < lists2.getLength(); i2++)
              {
                Node element2 =  lists2.item(i2);
                String type2 = element2.getNodeName();

                if (type2.equals("#text"))
                {
                  String value2=element2.getNodeValue();

                  parse_vars(value2,1);
                }
              }
            }
          }

          switch(modus)
          {
            case FORMAT_HTML:
              out.writeBytes("</table>\n");
              break;
            case FORMAT_TEX:
              out.writeBytes("\\end{supertabular}\n");
              out.writeBytes("\\end{center}\n\n");
              break;
          }

        }

        //===========================================================================

        if (type.equals("constants"))
        {
          switch(modus)
          {
            case FORMAT_HTML:
              out.writeBytes("\n<div class=\"constants\">\n    <a name=\"CONSTANT\">Constants</a><br>\n");
              out.writeBytes("<table border=\"0\">\n");
              out.writeBytes("    <th width=\"65\">Name</th>\n    <th>Type</th>\n    <th>Purpose</th>\n");
              break;
            case FORMAT_TEX:
              out.writeBytes("{\\bf Constant definitions: }\n");
              out.writeBytes("\\begin{center}\\small\n");
              out.writeBytes("\\phantom{fucksupertabular}\\begin{supertabular}{l l l p{5cm}}\n");
              out.writeBytes("Name&Type&&Purpose\\\\ \\hline \n");
              break;
          }

          NodeList lists1 = element.getChildNodes();
          for (int i1 = 0; i1 < lists1.getLength(); i1++)
          {
            Node element1 =  lists1.item(i1);
            String type1 = element1.getNodeName();

            if (type1.equals("constantblock"))
            {
              NamedNodeMap constantblockattr = element1.getAttributes();

              int varset = 0;
              String constdescr = "";
              String varname = "";

              for (int attridx = 0; attridx < constantblockattr.getLength(); attridx++)
              {
                if (constantblockattr.item(attridx).getNodeName().equals("variable"))
                {
                  varname=constantblockattr.item(attridx).getNodeValue();

                  varset=1;
                }

                if (constantblockattr.item(attridx).getNodeName().equals("description"))
                {
                  constdescr=constantblockattr.item(attridx).getNodeValue();
                }
              }

              if (varset==0)
              {
                switch(modus)
                {
                  case FORMAT_HTML:
                    out.writeBytes("<tr><td>&nbsp;</td></tr><tr><td colspan=3>Constants: </td></tr>\n");
                    if (!constdescr.equals("")) out.writeBytes("<tr><td>&nbsp;</td></tr><tr><td colspan=3>Purpose: "+constdescr+" </td></tr>\n");
                    break;
                  case FORMAT_TEX:
                    out.writeBytes("\\\\\n");
                    if (!constdescr.equals("")) out.writeBytes("Purpose: "+constdescr+"\\\\\n");
                    break;
                }
              }
              else
              {
                switch(modus)
                {
                  case FORMAT_HTML:
                    out.writeBytes("<tr><td>&nbsp;</td></tr><tr><td colspan=3><a name="+varname+"linkdef></a>Constants for: <a href=\"#"+varname+"\">"
                      +varname+"</a></td></tr>\n");
                    if (!constdescr.equals("")) out.writeBytes("<tr><td>&nbsp;</td></tr><tr><td colspan=3>Purpose: "+constdescr+" </td></tr>\n");
                    break;
                  case FORMAT_TEX:
                    out.writeBytes("Constants for "+varname+"\\\\\n");
                    if (!constdescr.equals("")) out.writeBytes("Purpose: "+constdescr+"\\\\\n");
                    break;
                }
              }

              NodeList lists2 = element1.getChildNodes();

              for (int i2 = 0; i2 < lists2.getLength(); i2++)
              {
                Node element2 =  lists2.item(i2);
                String type2 = element2.getNodeName();

                if (type2.equals("#text"))
                {
                  String value2=element2.getNodeValue();

                  parse_vars(value2,0);
                }
              }
            }
          }

          switch(modus)
          {
            case FORMAT_HTML:
              out.writeBytes("</table>\n");
              break;
            case FORMAT_TEX:
              out.writeBytes("\\end{supertabular}\n");
              out.writeBytes("\\end{center}\n\n");
              break;
          }
        }

        //===========================================================================

      }


      switch(modus)
      {
        case FORMAT_HTML:
          out.writeBytes("<br><br><hr><br><h2>Detailed subroutine and function documentation</h2>\n");
          break;
     //   case FORMAT_TEX:
     //     out.writeBytes("{\\bf Detailed subroutine and function documentation}\\\\[0.5cm]\n");
     //     break;
      }

      for (int i = 0; i < lists.getLength(); i++)
      {

        Node element =  lists.item(i);

        String type = element.getNodeName();

        if ((type.equals("subroutine")) || (type.equals("function")))
        {

          NodeList lists1 = element.getChildNodes();

          String result="";
          String collect;

          for (int i1 = 0; i1 < lists1.getLength(); i1++)
          {
            Node element1 =  lists1.item(i1);
            String type1 = element1.getNodeName();

            if (type1.equals("global"))
            {
              NamedNodeMap globalattr = element1.getAttributes();

              String varname = "";

              for (int attridx = 0; attridx < globalattr.getLength(); attridx++)
              {
                if (globalattr.item(attridx).getNodeName().equals("variable"))
                {
                  varname=globalattr.item(attridx).getNodeValue();
                }

              }

              switch(modus)
              {
                case FORMAT_HTML:
                  out.writeBytes("<br><b>Global variables:</b>"+varname+"<br>\n");
                  break;
                case FORMAT_TEX:
                  out.writeBytes("\\begin{sloppypar}\n{\\bf Global variables: }\\\\[0.1cm]\n");

                  String[] ents = varname.replaceAll("_","\\\\_").replaceAll(" ->"," \\\\ensuremath{\\\\rightarrow}").split(",");

                  for (int i2=0; i2<ents.length; i2++)
                    if (i2<ents.length-1)
                      out.writeBytes("{\\tt "+ents[i2]+",} ");
                    else
                      out.writeBytes("{\\tt "+ents[i2]+"} ");

                  out.writeBytes("\n\\end{sloppypar}\n\n");

                  break;
              }
            }


            if ((type1.equals("input"))||(type1.equals("output"))||(type1.equals("inoutput")))
            {
              NodeList lists2 = element1.getChildNodes();

              String value2=lists2.item(0).getNodeValue();

              if (type1.equals("input"))
              {
                switch(modus)
                {
                  case FORMAT_HTML:
                    out.writeBytes("<b>Input variables:</b><br>\n");
                    break;
                  case FORMAT_TEX:
                    out.writeBytes("\n{\\bf Input variables:}");
                    break;
                }
              }
              else
              {
                if (type1.equals("output"))
                {
                  switch(modus)
                  {
                    case FORMAT_HTML:
                      out.writeBytes("<br><b>Output variables:</b><br>\n");
                      break;
                    case FORMAT_TEX:
                      out.writeBytes("\n{\\bf Output variables:}");
                      break;
                   }
                }
                else
                {
                  switch(modus)
                  {
                    case FORMAT_HTML:
                      out.writeBytes("<br><b>Input/Output variables:</b><br>\n");
                      break;
                    case FORMAT_TEX:
                      out.writeBytes("\n{\\bf Input/Output variables:}");
                      break;
                   }
                }
              }

              switch(modus)
              {
                case FORMAT_HTML:
                  out.writeBytes("<table>\n");
                  break;
                case FORMAT_TEX:
                  out.writeBytes("\\begin{center}\\small\n");
                  out.writeBytes("\\phantom{fucksupertabular}\\begin{tabular}{l l l p{5cm}}\n");
                  out.writeBytes("Name&Type&Rank&Description \\\\ \\hline\n");
                  break;
              }

              parse_vars(value2,0);

              switch(modus)
              {
                case FORMAT_HTML:
                  out.writeBytes("</table>\n");
                  break;
                case FORMAT_TEX:
                  out.writeBytes("\\end{tabular}\n");
                  out.writeBytes("\\end{center}\n\n");
                  break;
              }
            }

            if (type1.equals("description"))
            {
              NodeList lists2 = element1.getChildNodes();

              String descr="";

              for (int i11 = 0; i11 < lists2.getLength(); i11++)
              {
                Node element2 =  lists2.item(i11);
                String type2 = element2.getNodeName();

                String textlines = "";

                if (type2.equals("#text"))
                {
                  if (modus==FORMAT_TEX)
                    textlines=lists2.item(i11).getNodeValue().replaceAll("_","\\\\_").replaceAll(" ->"," \\\\ensuremath{\\\\rightarrow}");
                  else
                    textlines=lists2.item(i11).getNodeValue();
                }
                else if (type2.equals("tex"))
                {
                  textlines="!"+lists2.item(i11).getChildNodes().item(0).getNodeValue();
                }
                else
                  System.out.println("Unknown description tag");


                String[] lines = textlines.split("\n");


                for(int i2=0;i2<lines.length;i2++)
                {
                  if (lines[i2].length()<=1)
                  {
                    if (i2>0)
                    {
                      switch(modus)
                      {
                        case FORMAT_HTML:
                          descr=descr.concat("<br><br>");
                          break;
              //          case FORMAT_TEX:
              //            descr=descr.concat("\n");
              //            break;
                      }
                    }
                  }
                  else
                  {
                    if(lines[i2].trim().length() > 0)
                    {
                      descr=descr.concat(lines[i2].trim().substring(1)+" ");
                    }
                    else
                    {
                      descr=descr.concat(" ");
                    }
                  }
                }

              }

              switch(modus)
              {
                case FORMAT_HTML:
                  out.writeBytes("<b>Description:</b><br>"+descr+"<br><br>\n");
                  break;
                case FORMAT_TEX:
                  out.writeBytes("\n{\\bf Description: }\\\\"+descr.replaceAll("%","\\\\%")+"\n");
                  break;
              }
            }



            if (type1.equals("result"))
            {
              NodeList lists2 = element1.getChildNodes();

              String[] lines = lists2.item(0).getNodeValue().split("\n");
              collect="";

              for(int l = 0; l<lines.length; l++)
              {
                if (lines[l].trim().indexOf("!")==0)
                {
                  collect=collect.concat(" "+lines[l].trim().substring(1));
                }
              }


              switch(modus)
              {
                case FORMAT_HTML:
                  out.writeBytes("<b>Result:</b><br> &nbsp;<em>"+result+"</em> &nbsp;"+collect.trim()+"<br><br>\n");
                  break;
                case FORMAT_TEX:
                    out.writeBytes("\n{\\bf Result:}\\\\[0.1cm] {\\it "
                    +result.replaceAll("_","\\\\_").replaceAll(" ->"," \\\\ensuremath{\\\\rightarrow}").replaceAll("recursive","")+"}: "
                    +collect.trim().replaceAll("_","\\\\_").replaceAll(" ->"," \\\\ensuremath{\\\\rightarrow}")+"\n");
                  break;
              }
            }

            if (type1.equals("#text"))
            {
              String value1=element1.getNodeValue();
	      if (patternSubroutine.matcher(value1).matches())
              {
                String mame=value1.substring(value1.toLowerCase().indexOf("subroutine")+10,value1.indexOf("(")).trim();

                switch(modus)
                {
                  case FORMAT_HTML:
                    out.writeBytes("<h3><a name="+mame+">Subroutine "
                     +value1.substring(value1.toLowerCase().indexOf("subroutine")+10).replaceAll("!","").trim()+"</a></h3>\n");
                    break;
                  case FORMAT_TEX:
                    String mame1=value1.substring(value1.toLowerCase().indexOf("subroutine")+10,value1.indexOf("("));
                    out.writeBytes("\\subsection{Subroutine "+mame1.replaceAll("_","\\\\_").replaceAll(" ->"," \\\\ensuremath{\\\\rightarrow}")+"}\n");

                    out.writeBytes("\\begin{sloppypar}\n{\\bf Interface: }\\\\[0.1cm]");

                    String iface = value1.substring(value1.toLowerCase().indexOf("subroutine")+10).replaceAll("!","").replaceAll("\n"," ").replaceAll("_","\\\\_").replaceAll(" ->"," \\\\ensuremath{\\\\rightarrow}").trim();

                   // System.out.println("TSCAK |"+iface+"|\n");

                    String[] ifaces = iface.split(",");

                    for(int l = 0; l<ifaces.length; l++)
                      if (l<ifaces.length-1)
                        out.writeBytes("{\\tt "+ifaces[l].trim()+"}, ");
                      else
                        out.writeBytes("{\\tt "+ifaces[l].trim()+"} ");

                    out.writeBytes("\n\\end{sloppypar}\n");

                    break;
                }
              }

	      if (patternFunction.matcher(value1).matches())
	      {
		String mame = value1.substring(value1.toLowerCase().indexOf("function") + 9, value1.lastIndexOf("(")).trim();

                result=value1.substring(0,value1.toLowerCase().indexOf("function"));

                switch(modus)
                {
                  case FORMAT_HTML:
                    out.writeBytes("<h3><a name="+mame+">Function "
                      +value1.substring(value1.toLowerCase().indexOf("function")+9).replaceAll("!","").trim()+"</a></h3>\n");
                    break;
                  case FORMAT_TEX:
                    String mame2=value1.substring(value1.toLowerCase().indexOf("function")+9);
                    String mame1=mame2.substring(0,mame2.indexOf("("));
                    out.writeBytes("\\subsection{Function "+mame1.replaceAll("_","\\\\_").replaceAll(" ->"," \\\\ensuremath{\\\\rightarrow}")+"}\n");
                    out.writeBytes("{\\bf Interface: }\\\\[0.1cm]{\\tt "+mame2.substring(0,mame2.indexOf(")")+1).replaceAll("!","").replaceAll("_","\\\\_").replaceAll(" ->"," \\\\ensuremath{\\\\rightarrow}")+"}\n");
                    break;
                }
              }
            }
          }

          switch(modus)
          {
            case FORMAT_HTML:
              out.writeBytes("<hr>\n");
              break;
          }
        }
      }


      switch(modus)
      {
        case FORMAT_HTML:

  //        out.writeBytes("</BODY>\n");
//        out.writeBytes("</HTML>\n");
          break;
      }

      out.close();
    }
    catch (Exception ex)
    {
      System.out.println("Exp "+ex);
      ex.printStackTrace();
          System.exit(1);
    }

  }
}
