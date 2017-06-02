// Copyright (c) 2012 The Foundry Visionmongers Ltd.  All Rights Reserved.

// Standard plug-in include files.
#include "DDImage/Iop.h"
#include "DDImage/Row.h"
#include "DDImage/Knobs.h"

using namespace DD::Image;

class ExamplePlugin : public Iop
{

public:

  //! Constructor. Initialize user controls to their default values.
  ExamplePlugin (Node* node) : Iop (node)
  {
  }

  ~ExamplePlugin () {}
  
  void _validate(bool);
  void _request(int x, int y, int r, int t, ChannelMask channels, int count);
  
  //! This function does all the work.
  void engine ( int y, int x, int r, ChannelMask channels, Row& outRow );

  //! Return the name of the class.
  const char* Class() const { return CLASS; }
  const char* node_help() const { return HELP; }

private:

  //! Information to the plug-in manager of DDNewImage/NUKE.
  static const Iop::Description description;
  static const char* const CLASS;
  static const char* const HELP;
}; 


/*! This is a function that creates an instance of the operator, and is
   needed for the Iop::Description to work.
 */
static Iop* ExamplePluginCreate(Node* node)
{
  return new ExamplePlugin(node);
}

/*! The Iop::Description is how NUKE knows what the name of the operator is,
   how to create one, and the menu item to show the user. The menu item may be
   0 if you do not want the operator to be visible.
 */
const Iop::Description ExamplePlugin::description ( CLASS, "Examples/ExamplePlugin",
                                                     ExamplePluginCreate );


const char* const ExamplePlugin::CLASS = "ExamplePlugin";
const char* const ExamplePlugin::HELP = "Example Plugin";

void ExamplePlugin::_validate(bool for_real)
{
  copy_info(0); // copy bbox channels etc from input0, which will validate it.
  info_.channels(Mask_RGB);

}

void ExamplePlugin::_request(int x, int y, int r, int t, ChannelMask channels, int count)
{
  // for this example, we're only interested in the RGB channels
  input(0)->request( x, y, r, t, ChannelMask(Mask_RGB), count );
}


/*! For each line in the area passed to request(), this will be called. It must
   calculate the image data for a region at vertical position y, and between
   horizontal positions x and r, and write it to the passed row
   structure. Usually this works by asking the input for data, and modifying
   it.

 */
void ExamplePlugin::engine ( int y, int x, int r,
                              ChannelMask channels, Row& outRow )
{
  // for this example, create a simple greyscale output from RGB channels
  ChannelMask rgbMask(Mask_RGB);

  Row inputRow(x, r);
  inputRow.get(input0(), y, x, r, rgbMask);

  for( int curX = x ; curX < r; curX++ ) {

    float newValue = 0.0f;

    foreach ( channel, rgbMask ) {
      newValue += inputRow[channel][curX];
    }
    newValue /= channels.size();

    foreach ( channel, channels ) {
      outRow.writable(channel)[curX] = newValue;
    }
  }
}
