

/*
This is a UI file (.ui.qml) that is intended to be edited in Qt Design Studio only.
It is supposed to be strictly declarative and only uses a subset of QML. If you edit
this file manually, you might introduce QML code that is not supported by Qt Design Studio.
Check out https://doc.qt.io/qtcreator/creator-quick-ui-forms.html for details on .ui.qml files.
*/
import QtQuick 2.15
import QtQuick.Controls 2.15
import QtQuick.Layouts

Item {
    ColumnLayout {
        id: columnLayout
        anchors.fill: parent
        spacing: 2

        Switch {
            id: switch1
            text: qsTr("Switch")
            onCheckedChanged: {
                if(!checked) resize_details.start()
                else resize_details_open.start()
            }
        }

        Column {
            id: column

            Layout.fillWidth: true

            Rectangle {
                id: rectangle
                color: "#00fa0000"
                clip: true
                height: column.height
                width: column.width

                SequentialAnimation {
                    id: resize_details
                    NumberAnimation {
                        target: rectangle;
                        property: "height";
                        to: 0;
                        duration: 100
                    }
                }
                SequentialAnimation {
                    id: resize_details_open
                    NumberAnimation {
                        target: rectangle;
                        property: "height";
                        to: 100;
                        duration: 100
                    }
                }

                GridLayout {
                    id: gridLayout
                    anchors.fill: parent
                    anchors.leftMargin: 5
                    anchors.rightMargin: 5
                    anchors.topMargin: 0
                    anchors.bottomMargin: 0
                    rows: 2
                    columns: 3

                    Label {
                        id: label1
                        text: qsTr("Label")
                        horizontalAlignment: Text.AlignHCenter
                        Layout.fillWidth: true
                    }

                    Slider {
                        id: slider
                        value: 0.5
                        Layout.fillWidth: true
                    }


                    Label {
                        id: label3
                        text: qsTr("Label")
                        horizontalAlignment: Text.AlignHCenter
                        Layout.fillWidth: true
                    }

                    Label {
                        id: label
                        text: qsTr("Label")
                        horizontalAlignment: Text.AlignHCenter
                        Layout.fillWidth: true
                    }


                    Slider {
                        id: slider1
                        value: 0.5
                        Layout.fillWidth: true
                    }


                    Label {
                        id: label2
                        text: qsTr("Label")
                        horizontalAlignment: Text.AlignHCenter
                        Layout.fillWidth: true
                    }

                    Item {
                        id: item2
                        Layout.fillHeight: true
                        Layout.fillWidth: true
                    }



                }
            }
        }

        Rectangle {
            id: rectangle1
            width: 200
            height: 200
            color: "#0ae14e"
            Layout.fillHeight: true
            Layout.fillWidth: true
        }

        Item {
            id: item1
            Layout.fillHeight: true
            Layout.fillWidth: true
        }

    }
}
